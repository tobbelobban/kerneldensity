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

#include <KTH/kerneldensityestimation/processors/sequentialkde.h>

#include <inviwo/core/util/glm.h>
#include <inviwo/core/datastructures/volume/volume.h>
#include <inviwo/core/datastructures/volume/volumeramprecision.h>

#include <inviwo/core/datastructures/buffer/buffer.h>
#include <inviwo/core/datastructures/buffer/bufferramprecision.h>
#include <inviwo/core/datastructures/buffer/bufferram.h>

#include <omp.h>

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo SequentialKDE::processorInfo_{
    "org.inviwo.SequentialKDE",      // Class identifier
    "Sequential KDE",                // Display name
    "Undefined",              // Category
    CodeState::Experimental,  // Code state
    Tags::None,               // Tags
};
const ProcessorInfo SequentialKDE::getProcessorInfo() const { return processorInfo_; }

SequentialKDE::SequentialKDE()
    : Processor()
    , volume_in_("volumeInport")
	, mesh_in_("meshInport")
	, volume_out_("volumeOutport")
	, bandwidth_prop("KDEbandwidth", "Bandwidth", 0.0001, 0.00001, 3, 0.00001)
	, cutoff_prop("cutoff", "Cutoff", 0.0001, 0.00001, 3, 0.00001)
	, fast_KDE_prop("fast_kde", "Use fast KDE", true)    
{
    
	addPort(volume_in_);
    addPort(mesh_in_);
    addPort(volume_out_);

	addProperty(bandwidth_prop);
	addProperty(fast_KDE_prop);
	addProperty(cutoff_prop);
    
}


void SequentialKDE::KDE() {
    auto buff_pair = mesh_in_.getData()->findBuffer(BufferType::IndexAttrib);
    const size3_t vol_dims = volume_in_.getData()->getDimensions();

    if (buff_pair.first == nullptr) {
        kde_vol = std::make_shared<Volume>(vol_dims, DataFloat32::get());
        volume_out_.setData(kde_vol);
		LogProcessorError("ERROR: No index buffer found.");
        return;
    }

	kde_vol = std::make_shared<Volume>(vol_dims, DataFloat32::get());
	double maxval = 0;

	// compute & store KDE
    kde_vol->getEditableRepresentation<VolumeRAM>()->dispatch<void, dispatching::filter::Scalars>(
		[&](auto kdevol_pr) {
			const size_t vol_size = vol_dims.x * vol_dims.y * vol_dims.z;
			const auto basis = volume_in_.getData()->getBasis();
			const dvec3 offset = volume_in_.getData()->getOffset();
			const dvec3 spacing(basis[0][0]/(vol_dims.x-1), basis[1][1]/(vol_dims.y-1), basis[2][2]/(vol_dims.z-1));
			LogProcessorInfo(spacing);
			const size_t extrema_size = buff_pair.first->getRepresentation<BufferRAM>()->getSize();
			const int* extrema_idx = static_cast<const int*>(buff_pair.first->getRepresentation<BufferRAM>()->getData());
	
			const double h = bandwidth_prop.get();
			
			const int xy_ = vol_dims.x * vol_dims.y;
			const double denom = (h*h*h*(double)extrema_size);
	
			using KDEType = util::PrecisionValueType<decltype(kdevol_pr)>;
			KDEType* kde_data = kdevol_pr->getDataTyped();

			for (int i = 0; i < extrema_size; ++i) {
				const dvec3 extrema_pos = offset + spacing * dvec3(extrema_idx[i] % vol_dims.x, (extrema_idx[i] / vol_dims.x) % vol_dims.y, extrema_idx[i] / xy_);
				#pragma omp parallel for shared(extrema_pos, kde_data, spacing, offset, h)
				for (int iz = 0; iz < vol_dims.z; ++iz) {
					int z_offset = iz * xy_;
					for (int iy = 0; iy < vol_dims.y; ++iy) {
						int arr_offset = iy * vol_dims.x + z_offset;
						for (int ix = 0; ix < vol_dims.x; ++ix) {
							dvec3 curr_world_pos = dvec3(ix, iy, iz) * spacing + offset;
							dvec3 offset_from_extrema = (curr_world_pos - extrema_pos)/h;
							double r_sq = offset_from_extrema.x * offset_from_extrema.x +
										  offset_from_extrema.y * offset_from_extrema.y +
										  offset_from_extrema.z * offset_from_extrema.z;
							kde_data[arr_offset++] += (1.0 / (2.0 * M_PI * sqrt(2.0 * M_PI))) * exp(-0.5 * r_sq);
						}
					}
				}
			}
			int index = 0;
			for (int iz = 0; iz < vol_dims.z; ++iz) {
				for (int iy = 0; iy < vol_dims.y; ++iy) {
					for (int ix = 0; ix < vol_dims.x; ++ix) {
						kde_data[index] /= denom;
						if(kde_data[index] > maxval) maxval = kde_data[index];
						++index;
					}
				}
			}
		}
	);	

	kde_vol->copyMetaDataFrom(*volume_in_.getData());
	kde_vol->dataMap_ = volume_in_.getData()->dataMap_;
    kde_vol->dataMap_.dataRange = vec2(0,maxval);
	kde_vol->dataMap_.valueRange = vec2(0,maxval);
	
	kde_vol->setBasis(volume_in_.getData()->getBasis());
    kde_vol->setOffset(volume_in_.getData()->getOffset());
    kde_vol->setWorldMatrix(volume_in_.getData()->getWorldMatrix());
    kde_vol->setModelMatrix(volume_in_.getData()->getModelMatrix());
	
	volume_out_.setData(kde_vol);
}

void SequentialKDE::makeKDEStencil(const int extrema_size, const double h) {
	const double cutoff = cutoff_prop.get();
	const size3_t vol_dims = volume_in_.getData()->getDimensions();
	const auto basis = volume_in_.getData()->getBasis();
	
	const dvec3 spacing(basis[0][0]/(vol_dims.x-1), basis[1][1]/(vol_dims.y-1), basis[2][2]/(vol_dims.z-1));

	stencil_half_dims = std::make_shared<ivec3>(std::ceil(cutoff / spacing.x),
												std::ceil(cutoff / spacing.y),
												std::ceil(cutoff / spacing.z));

    stencil_dims = std::make_shared<ivec3>(	stencil_half_dims->x * 2 + 1, 
											stencil_half_dims->y * 2 + 1,
											stencil_half_dims->z * 2 + 1 );
	
	const int stencil_size = stencil_dims->x * stencil_dims->y * stencil_dims->z;
	kde_stencil = std::make_shared<std::vector<float>>();
	kde_stencil->reserve(stencil_size);	
	
	double val;
	const double denom = (h*h*h * extrema_size);
	const dvec3 scaled_spacing = spacing / h;
	int count = 0;
	double sum = 0;
	for (int iz = -stencil_half_dims->z; iz < stencil_half_dims->z + 1; ++iz) {
		for (int iy = -stencil_half_dims->y; iy < stencil_half_dims->y + 1; ++iy) {
			for (int ix = -stencil_half_dims->x; ix < stencil_half_dims->x + 1; ++ix) {				
                double r_sq =	pow( ix * scaled_spacing.x, 2) +
								pow( iy * scaled_spacing.y, 2) +
								pow( iz * scaled_spacing.z, 2);
				val = (1.0 / (2.0 * M_PI * sqrt(2.0 * M_PI))) * exp(-0.5 * r_sq);
				sum += val;
                kde_stencil->push_back( val/denom );
				count += 1;
			}
		}
	}
	double mean = sum / count;
	LogProcessorInfo("Mean:");
	LogProcessorInfo(mean);
	LogProcessorInfo("Sum:");
	LogProcessorInfo(sum);
	LogProcessorInfo("Denom:");
	LogProcessorInfo(denom);
}

void SequentialKDE::fastKDE() {
    auto buff_pair = mesh_in_.getData()->findBuffer(BufferType::IndexAttrib);
    const size3_t vol_dims = volume_in_.getData()->getDimensions();
    const size_t vol_size = vol_dims.x * vol_dims.y * vol_dims.z;

    if (buff_pair.first == nullptr) {
		LogProcessorInfo("Error: no index buffer.");
        kde_vol = std::make_shared<Volume>(vol_dims, DataFloat32::get());
        volume_out_.setData(kde_vol);
        return;
    }

    const size_t extrema_size = buff_pair.first->getRepresentation<BufferRAM>()->getSize();
	LogProcessorInfo("Number extrema:");
	LogProcessorInfo(extrema_size);
    const int* extrema_idx = static_cast<const int*>(buff_pair.first->getRepresentation<BufferRAM>()->getData());
    const double h = bandwidth_prop.get();
	
	// make stencil for fast KDE
	makeKDEStencil(extrema_size, h);
	kde_vol = std::make_shared<Volume>(vol_dims, DataFloat32::get());
    volume_out_.setData(kde_vol);
    
	// compute & store KDE
    kde_vol = std::make_shared<Volume>(vol_dims, DataFloat32::get());
    double maxval = 0;

    kde_vol->getEditableRepresentation<VolumeRAM>()->dispatch<void, dispatching::filter::Scalars>(
		[&](auto kdevol_pr) {
        
			using KDEType = util::PrecisionValueType<decltype(kdevol_pr)>;
			KDEType* kde_data = kdevol_pr->getDataTyped();
		
			const int xy_ = vol_dims.x * vol_dims.y;
			const int stencil_xy_ = stencil_dims->x * stencil_dims->y;			
			double * max_values = new double[omp_get_max_threads()]{0};

 			for (int i = 0; i < extrema_size; ++i) {
				const size3_t extrema_idx(extrema_idx[i] % vol_dims.x, (extrema_idx[i] / vol_dims.x) % vol_dims.y, extrema_idx[i] / xy_);
				#pragma omp parallel for default(shared)
				for (int iz = extrema_idx.z - stencil_half_dims->z; iz <= extrema_idx.z + stencil_half_dims->z; ++iz) {
					if (iz < 0) continue;
					if (iz >= vol_dims.z) break;
					const size_t stencil_iz = iz - extrema_idx.z + stencil_half_dims->z;
					for (int iy = extrema_idx.y - stencil_half_dims->y; iy <= extrema_idx.y + stencil_half_dims->y; ++iy) {
						if (iy < 0) continue;	
						if(iy >= vol_dims.y) break;
						const size_t stencil_iy = iy - extrema_idx.y + stencil_half_dims->y;
						for (int ix = extrema_idx.x - stencil_half_dims->x; ix <= extrema_idx.x + stencil_half_dims->x; ++ix) {
							if (ix < 0) continue;
							if(ix >= vol_dims.x) break;
							const size_t curr_i = ix + iy * vol_dims.x + iz * xy_;
							const size_t stencil_i = ix - extrema_idx.x + stencil_half_dims->x + stencil_iy * stencil_dims->x + stencil_iz * stencil_xy_;
							kde_data[curr_i] += KDEType(kde_stencil->at(stencil_i));
							if(i < extrema_size-1) continue; // wait until last extrema to compute maximum
							if(kde_data[curr_i] > max_values[omp_get_thread_num()]) max_values[omp_get_thread_num()] = kde_data[curr_i];
						}
					}
				}
			}
			for(int i = 0; i < omp_get_max_threads(); ++i) {
				if(max_values[i] > maxval) maxval = max_values[i];
			}
			delete[] max_values;
		}
	);

	kde_vol->copyMetaDataFrom(*volume_in_.getData());
	kde_vol->dataMap_ = volume_in_.getData()->dataMap_;
    kde_vol->dataMap_.dataRange = vec2(0,maxval);
	kde_vol->dataMap_.valueRange = vec2(0,maxval);
	
	kde_vol->setBasis(volume_in_.getData()->getBasis());
    kde_vol->setOffset(volume_in_.getData()->getOffset());
    kde_vol->setWorldMatrix(volume_in_.getData()->getWorldMatrix());
    kde_vol->setModelMatrix(volume_in_.getData()->getModelMatrix());
	
	volume_out_.setData(kde_vol);
}

void SequentialKDE::process() {
	
	auto start = std::chrono::high_resolution_clock::now();
	
	if(fast_KDE_prop.get()) {
		fastKDE();
	} else {
		KDE();
	}
    
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::string time_str = "KDE: " + std::to_string(duration.count()) + " ms";
    
	LogProcessorInfo(time_str);
}

}  // namespace inviwo
