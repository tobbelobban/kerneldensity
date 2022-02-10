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

#include <inviwo/kerneldensityestimation/processors/slicesequentialkde.h>

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo SliceSequentialKDE::processorInfo_{
    "org.inviwo.SliceSequentialKDE",      // Class identifier
    "Slice Sequential KDE",                // Display name
    "Undefined",              // Category
    CodeState::Experimental,  // Code state
    Tags::None,               // Tags
};
const ProcessorInfo SliceSequentialKDE::getProcessorInfo() const { return processorInfo_; }

SliceSequentialKDE::SliceSequentialKDE()
    : Processor()
    , volume_inport_("volume_inport_")
	, mesh_inport_("mesh_inport_")
	, volume_outport_("volume_outport_")
	, bandwidth_prop_("KDEbandwidth", "Bandwidth", 0.0001f, 0.00001f, 3.0f, 0.00001f)
	, dimension_("dimension", "Dimension",
              {{"x_dim", "X", Dimension::X},
               {"y_dim", "Y", Dimension::Y},
               {"z_dim", "Z", Dimension::Z}}) {

	addPort(volume_inport_);
	addPort(mesh_inport_);
    addPort(volume_outport_);
    
	addProperty(dimension_);
	addProperty(bandwidth_prop_);
}

void SliceSequentialKDE::process() {
	// for timing processor
	auto start = std::chrono::high_resolution_clock::now();
    const std::shared_ptr<const Volume> volume_in_ptr = volume_inport_.getData();
	auto buff_pair = mesh_inport_.getData()->findBuffer(BufferType::IndexAttrib);
	const size3_t vol_dims = volume_in_ptr->getDimensions();
	const auto vol_xy = vol_dims.x * vol_dims.y;
	const ivec3 delta = {1, vol_dims.x, vol_xy};
	// compute grid spacing
	const vec3 spacing_offset = volume_in_ptr->getOffset();
	const auto basis = volume_in_ptr->getBasis();
	const vec3 spacing(	basis[0][0]/(vol_dims.x-1),
						basis[1][1]/(vol_dims.y-1),
						basis[2][2]/(vol_dims.z-1)	);
	// get the dimension to slice along 
	const Dimension slice_dim = dimension_.get();
	const Dimension dim1 = (slice_dim == Dimension::X ? Dimension::Y : Dimension::X);
	const Dimension dim2 = (dim1 == Dimension::X ? (slice_dim == Dimension::Y ? Dimension::Z : Dimension::Y) : Dimension::Z);
	const int islice_dim = static_cast<const int>(slice_dim);
	const int idim1 = static_cast<const int>(dim1);
	const int idim2 = static_cast<const int>(dim2);
	// set up write-to volume
	auto sliced_KDE_vol_repr = std::make_shared<VolumeRAMPrecision<float>>(vol_dims);
    float* sliced_KDE_raw_ptr = sliced_KDE_vol_repr->getDataTyped();
	// get extrema indices
	const size_t num_extrema = buff_pair.first->getRepresentation<BufferRAM>()->getSize();
	const int* extrema_indices = static_cast<const int*>(buff_pair.first->getRepresentation<BufferRAM>()->getData());
	// precompute constants
	const float KDE_constant = static_cast<float>(1.0f / (2.0f * M_PI));
	const double h = bandwidth_prop_.get();
	// counters for number of extrema per slice
	int* const extrema_per_slice_counter = new int[vol_dims[islice_dim]]{0};
	//compute sliced KDE
	for(int i = 0; i < num_extrema; ++i) {
		// 3D index coords for i'th extreme value 
		const ivec3 ev_grid_coord(	extrema_indices[i] % vol_dims.x,				// x
									(extrema_indices[i] / vol_dims.x) % vol_dims.y,	// y
									extrema_indices[i] / vol_xy	);					// z
		// 2D physical pos of current extreme value
		const vec2 extrema_pos = vec2(spacing[idim1], spacing[idim2]) * vec2(ev_grid_coord[idim1], ev_grid_coord[idim2]);
		extrema_per_slice_counter[ev_grid_coord[islice_dim]]++;
		// offset to start of slice in flattened array
		int offset = ev_grid_coord[islice_dim] * delta[islice_dim];
		// iterate over slice
		for(int i2 = 0; i2 < vol_dims[idim2]; ++i2) {
			for(int i1 = 0; i1 < vol_dims[idim1]; ++i1) {
				vec2 curr_world_pos = vec2(i1,i2) * vec2(spacing[idim1], spacing[idim2]);
				vec2 offset_from_ev = (curr_world_pos - extrema_pos)/h;
				double r_sq =	static_cast<double>(offset_from_ev.x) * offset_from_ev.x +
								offset_from_ev.y * offset_from_ev.y;
				sliced_KDE_raw_ptr[offset + i1*delta[idim1]] += static_cast<float>(KDE_constant * exp(-0.5 * r_sq));
			}	
			offset += delta[idim2];
		}
	}
	// scale by number of extrema
	const double denom = h*h;
	int index = 0;
	float max_val = 0.0;
	for (int iz = 0; iz < vol_dims.z; ++iz) {
		for (int iy = 0; iy < vol_dims.y; ++iy) {
			for (int ix = 0; ix < vol_dims.x; ++ix) {
				sliced_KDE_raw_ptr[index] /= static_cast<float>((extrema_per_slice_counter[ivec3(ix,iy,iz)[islice_dim]] * denom));
				if(sliced_KDE_raw_ptr[index] > max_val) max_val = sliced_KDE_raw_ptr[index];
				++index;
			}
		}
	}
	// free mem allocation
	delete[] extrema_per_slice_counter;
	// create output volume
	std::shared_ptr<Volume> KDE_sliced_vol = std::make_shared<Volume>(sliced_KDE_vol_repr);
    KDE_sliced_vol->setBasis(volume_in_ptr->getBasis());
	KDE_sliced_vol->setOffset(volume_in_ptr->getOffset());
	KDE_sliced_vol->copyMetaDataFrom(*volume_in_ptr);	
	KDE_sliced_vol->dataMap_.valueRange = vec2(0, max_val);
	KDE_sliced_vol->dataMap_.dataRange = vec2(0, max_val);
	KDE_sliced_vol->setModelMatrix(volume_in_ptr->getModelMatrix());
	KDE_sliced_vol->setWorldMatrix(volume_in_ptr->getWorldMatrix());
	//output time
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::string time_str = "KDE: " + std::to_string(duration.count()) + " ms";
	LogProcessorInfo(time_str);
	// set data on outport
    volume_outport_.setData(KDE_sliced_vol);
}

}  // namespace inviwo
