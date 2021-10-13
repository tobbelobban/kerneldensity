/*********************************************************************************
 *
 * Inviwo - Interactive Visualization Workshop
 *
 * Copyright (c) 2021 Inviwo Foundation
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *********************************************************************************/

#include <KTH/kerneldensityestimation/processors/fullkdecl.h>

#include <inviwo/core/util/glm.h>
#include <inviwo/core/datastructures/volume/volume.h>
#include <inviwo/core/datastructures/volume/volumeramprecision.h>

#include <inviwo/core/datastructures/buffer/buffer.h>
#include <inviwo/core/datastructures/buffer/bufferramprecision.h>
#include <inviwo/core/datastructures/buffer/bufferram.h>

#include <modules/opencl/kernelowner.h>
#include <modules/opencl/syncclgl.h>
#include <modules/opencl/buffer/buffercl.h>
#include <modules/opencl/buffer/bufferclconverter.h>
#include <modules/opencl/volume/volumecl.h>

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo FullKDECL::processorInfo_{
    "org.inviwo.FullKDECL",      // Class identifier
    "Full KDECL",                // Display name
    "Undefined",              // Category
    CodeState::Experimental,  // Code state
    Tags::None,               // Tags
};
const ProcessorInfo FullKDECL::getProcessorInfo() const { return processorInfo_; }

FullKDECL::FullKDECL()
    : Processor()
    , ProcessorKernelOwner(this)
	, volume_in_("volume_in")
	, mesh_in_("mesh_in")
    , volume_out_("volume_out")
	, bandwidth_prop("KDEbandwidth", "Bandwidth", 0.0001, 0.00001, 3, 0.00001)
	, cutoff_prop("cutoff", "Cutoff", 0.0001, 0.00001, 3, 0.00001) {

    addPort(volume_in_);
	addPort(mesh_in_);
    addPort(volume_out_);

	addProperty(bandwidth_prop);
	addProperty(cutoff_prop);

    KDEkernel_ = addKernel("kdekernel.cl", "volumeKDEKernel", "", "");
    
}

void FullKDECL::makeKDEStencilBuffer(const int nr_extrema, const double h, const double cutoff) {
	const size3_t vol_dims = volume_in_.getData()->getDimensions();
	const auto basis = volume_in_.getData()->getBasis();
	const dvec3 spacing(basis[0][0]/(vol_dims.x-1), basis[1][1]/(vol_dims.y-1), basis[2][2]/(vol_dims.z-1));

	stencil_half_dims = std::make_shared<ivec3>(std::ceil(cutoff / spacing.x),
												std::ceil(cutoff / spacing.y),
												std::ceil(cutoff / spacing.z));
    stencil_dims = std::make_shared<ivec3>(	stencil_half_dims->x * 2 + 1, 
											stencil_half_dims->y * 2 + 1,
											stencil_half_dims->z * 2 + 1);
	
	const int stencil_size = stencil_dims->x * stencil_dims->y * stencil_dims->z;

	stencil_buffer = std::shared_ptr<Buffer<float>>(new Buffer<float>(stencil_size));

	stencil_buffer->getEditableRepresentation<BufferRAM>()->dispatch<void, dispatching::filter::Scalars>(
		[&](auto stencil_pr) {
			using StencilType = util::PrecisionValueType<decltype(stencil_pr)>;
			std::vector<StencilType>& stencil_data = stencil_pr->getDataContainer();
			double val;
			int index = 0;
			for (int iz = -stencil_half_dims->z; iz < stencil_half_dims->z + 1; ++iz) {
				for (int iy = -stencil_half_dims->y; iy < stencil_half_dims->y + 1; ++iy) {
					for (int ix = -stencil_half_dims->x; ix < stencil_half_dims->x + 1; ++ix) {
						dvec3 curr_pos(ix * spacing.x / h, iy * spacing.y / h, iz * spacing.z / h);						
						double r_sq = curr_pos.x * curr_pos.x +
									curr_pos.y * curr_pos.y +
									curr_pos.z * curr_pos.z;
				
						val = (1.0 / (2.0 * M_PI * sqrt(2.0 * M_PI))) * exp(-0.5 * r_sq);
						stencil_data[index++] = StencilType(val/(h * nr_extrema));
					}
				}
			}
		}
	);
}
double FullKDECL::findVolumeMax() {
	const auto vol_dims = volume_in_.getData()->getDimensions();
    const int vol_size = vol_dims.x * vol_dims.y * vol_dims.z;
	return volume_out_.getData()->getRepresentation<VolumeRAM>()->dispatch<double, dispatching::filter::Scalars>(
		[&](auto out_vol_pr) {
			using VolType = util::PrecisionValueType<decltype(out_vol_pr)>;
			const VolType* vol_data = out_vol_pr->getDataTyped();
			VolType max = VolType(0);
			int index = 0;
			for (int i = 0; i < vol_size; ++i) {
					if(vol_data[i] > max) max = vol_data[i];
			}
			return max;
		}	
	);
}

void FullKDECL::process() {
	auto buff_pair = mesh_in_.getData()->findBuffer(BufferType::IndexAttrib);

	if (buff_pair.first == nullptr) {
		LogProcessorInfo("No extrema buffer");
        return;
    }

	const int* extrema = static_cast<const int*>(buff_pair.first->getRepresentation<BufferRAM>()->getData());
	const size_t extrema_size = buff_pair.first->getRepresentation<BufferRAM>()->getSize();    

	makeKDEStencilBuffer(extrema_size, bandwidth_prop.get(), cutoff_prop.get());

	const auto vol_dims = volume_in_.getData()->getDimensions();
    const int vol_size = vol_dims.x * vol_dims.y * vol_dims.z;

	if (!KDEkernel_) { // failed to load OpenCL kernel
		 volume_out_.setData(std::make_shared<Volume>(vol_dims, DataFloat32::get()));
		 LogProcessorInfo("Kernel build error!");
		 return;
	}
	
	auto out_vol = std::make_shared<Volume>(vol_dims, DataFloat32::get());
    out_vol->setModelMatrix(volume_in_.getData()->getModelMatrix());
    out_vol->setWorldMatrix(volume_in_.getData()->getWorldMatrix());
	out_vol->copyMetaDataFrom(*volume_in_.getData());
	out_vol->dataMap_ = volume_in_.getData()->dataMap_;
	
	out_vol->setBasis(volume_in_.getData()->getBasis());
    out_vol->setOffset(volume_in_.getData()->getOffset());

    volume_out_.setData(out_vol);

    VolumeCL* volumeOutCL = out_vol->getEditableRepresentation<VolumeCL>();
	auto tmpVolume_ = std::unique_ptr<Buffer<float>>(new Buffer<float>(vol_size));
	auto tmpVolumeCL = tmpVolume_->getEditableRepresentation<BufferCL>();
	auto stencilCL = stencil_buffer->getRepresentation<BufferCL>();
	auto extremaCL = buff_pair.first->getRepresentation<BufferCL>();
	
	LogProcessorInfo("Stencil dims:");
	LogProcessorInfo(*stencil_dims);
	const int xy_ = vol_dims.x * vol_dims.y;
	for(int i = 0; i < extrema_size; ++i) {
		cl::Event events[2];
		ivec3 extrema_idx(extrema[i] % vol_dims.x, (int)(extrema[i] / vol_dims.x) % vol_dims.y, extrema[i] / xy_);
		try {
			int kernel_idx = 0;
			KDEkernel_->setArg(kernel_idx++, *stencilCL);
			KDEkernel_->setArg(kernel_idx++, *stencil_dims);
			KDEkernel_->setArg(kernel_idx++, *stencil_half_dims);
			KDEkernel_->setArg(kernel_idx++, ivec3(vol_dims));
			//KDEkernel_->setArg(kernel_idx++, *extremaCL);
			KDEkernel_->setArg(kernel_idx++, extrema_idx);
			KDEkernel_->setArg(kernel_idx++, *tmpVolumeCL);
			auto res1 = OpenCL::getPtr()->getQueue().enqueueNDRangeKernel(
				*KDEkernel_, 
				cl::NullRange,  
				size3_t(*stencil_dims), 
				size3_t(1,1,1), 
				nullptr, 
				&events[0] );	
			std::vector<cl::Event> waitFor(1, events[0]);
			if(i != extrema_size-1) {
				cl::WaitForEvents(waitFor); // wait for each iteration to complete before queuing another
			} else {
				auto res2 = OpenCL::getPtr()->getQueue().enqueueCopyBufferToImage(
					tmpVolumeCL->get(), 
					volumeOutCL->getEditable(), 
					0,
					size3_t(0),
					size3_t(vol_dims),
					&waitFor, &events[1] );
				std::vector<cl::Event> waitFor2(1, events[1]);
				cl::WaitForEvents(waitFor2);
			}
		} catch (cl::Error& err) {
			LogError(getCLErrorString(err));
		}
	}

	stencil_buffer->removeRepresentation(stencilCL);
	auto max_val = findVolumeMax();
	out_vol->dataMap_.dataRange = vec2(0,max_val);
	out_vol->dataMap_.valueRange = vec2(0,max_val);
}

}  // namespace inviwo