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
	, use_scaling_prop("scaling", "Use scaling", true)  
{
    
	addPort(volume_in_);
    addPort(mesh_in_);
    addPort(volume_out_);

	addProperty(bandwidth_prop);
	addProperty(fast_KDE_prop);
	addProperty(cutoff_prop);
	addProperty(use_scaling_prop);
    
}


void SequentialKDE::KDE() {
	// det volume properties
    auto buff_pair = mesh_in_.getData()->findBuffer(BufferType::IndexAttrib);
    const size3_t vol_dims = volume_in_.getData()->getDimensions();
	const size_t vol_size = vol_dims.x * vol_dims.y * vol_dims.z;
	const auto basis = volume_in_.getData()->getBasis();
	const dvec3 spacing(basis[0][0]/(vol_dims.x-1), basis[1][1]/(vol_dims.y-1), basis[2][2]/(vol_dims.z-1));
	// check extreme value buffer
    if (buff_pair.first == nullptr) {
        LogProcessorError("ERROR: No index buffer found.");
		volume_out_.setData(std::make_shared<Volume>(vol_dims, DataFloat32::get()));
        return;
    }
	const size_t extrema_size = buff_pair.first->getRepresentation<BufferRAM>()->getSize();
	const int* extrema_idx = static_cast<const int*>(buff_pair.first->getRepresentation<BufferRAM>()->getData());
	// init output volume
	auto KDE_vr_precision = std::make_shared<VolumeRAMPrecision<float>>(vol_dims);
    float* KDE_raw_ptr = KDE_vr_precision->getDataTyped();
	// for finding maximum value of KDE
	double max_val = 0;
	// get bandwidth
	const double h = bandwidth_prop.get();
	// Explicitly disable dynamic teams
	omp_set_dynamic(0);     
	// Use 8 threads for all consecutive parallel regions
	omp_set_num_threads(16); 
	// compute & store KDE
    const int xy_ = vol_dims.x * vol_dims.y;
	const long double KDE_constant = (2.0 * M_PI * sqrt(2.0 * M_PI));
	const long double denom = (pow(h,3) * extrema_size);
	for (int i = 0; i < extrema_size; ++i) {
		const dvec3 extrema_pos = spacing * dvec3(	extrema_idx[i] % vol_dims.x,
													(extrema_idx[i] / vol_dims.x) % vol_dims.y,
													extrema_idx[i] / xy_	);
		// each thread iterates over a slice of the volume and applies KDE
		#pragma omp parallel shared(extrema_pos, KDE_raw_ptr, spacing, vol_dims, h, xy_)
		{
			double local_max = max_val;
			#pragma omp for nowait
			for (int iz = 0; iz < vol_dims.z; ++iz) {
				int buff_offset = iz * xy_;
				for (int iy = 0; iy < vol_dims.y; ++iy) {
					for (int ix = 0; ix < vol_dims.x; ++ix) {
						// calculate current physical pos of grid vertex
						dvec3 curr_world_pos = dvec3(ix, iy, iz) * spacing;
						// get distance from current extremum
						dvec3 scaled_delta_pos = (curr_world_pos - extrema_pos)/h;
						// apply gaussian kernel scaled by h
						double r_sq =	(scaled_delta_pos.x * scaled_delta_pos.x +
										scaled_delta_pos.y * scaled_delta_pos.y +
										scaled_delta_pos.z * scaled_delta_pos.z);
						KDE_raw_ptr[buff_offset] += KDE_constant * exp(-0.5 * r_sq); // add weight from Gaussian kernel
						if(KDE_raw_ptr[buff_offset] > local_max) local_max = KDE_raw_ptr[buff_offset];
					}
				}
			}
			#pragma omp critical
			{
				if(local_max > max_val) {
					max_val = local_max;
				}
			}
		}
	}
	// iterate over entire volume and divide by denom
	if(use_scaling_prop.get()) {
	int index = 0;
	for (int iz = 0; iz < vol_dims.z; ++iz) {
		for (int iy = 0; iy < vol_dims.y; ++iy) {
			for (int ix = 0; ix < vol_dims.x; ++ix) {
				KDE_raw_ptr[index] /= denom;
				++index;
			}
		}
	}
	max_val /= denom;
	}
	
	auto kde_vol = std::make_shared<Volume>(KDE_vr_precision);
	kde_vol->copyMetaDataFrom(*volume_in_.getData());
	kde_vol->dataMap_ = volume_in_.getData()->dataMap_;
    kde_vol->dataMap_.dataRange = vec2(0,max_val);
	kde_vol->dataMap_.valueRange = vec2(0,max_val);
	kde_vol->setBasis(volume_in_.getData()->getBasis());
    kde_vol->setOffset(volume_in_.getData()->getOffset());
    kde_vol->setWorldMatrix(volume_in_.getData()->getWorldMatrix());
    kde_vol->setModelMatrix(volume_in_.getData()->getModelMatrix());
	volume_out_.setData(kde_vol);
}

void SequentialKDE::makeKDEStencil(const int extrema_size, const double h) {
	const double cutoff = cutoff_prop.get();
	const ivec3 vol_dims = volume_in_.getData()->getDimensions();
	const auto basis = volume_in_.getData()->getBasis();
	// compute grid spacing
	const dvec3 spacing(	basis[0][0]/(vol_dims.x-1),
							basis[1][1]/(vol_dims.y-1),
							basis[2][2]/(vol_dims.z-1)	);
	// store half the stencil size
	stencil_half_dims = std::make_shared<ivec3>(	std::ceil(cutoff / spacing.x),
													std::ceil(cutoff / spacing.y),
													std::ceil(cutoff / spacing.z)	);
	// store the entire stencil size
    stencil_dims = std::make_shared<ivec3>(	stencil_half_dims->x * 2 + 1, 
											stencil_half_dims->y * 2 + 1,
											stencil_half_dims->z * 2 + 1	);
	// ensure stencil isnt larger than volume dimensions
	/*if(stencil_dims->x > vol_dims.x) {
		if(vol_dims.x % 2) {
			stencil_half_dims->x = vol_dims.x / 2;
			stencil_dims->x = vol_dims.x;
		} else {
			stencil_half_dims->x = vol_dims.x / 2 - 1;
			stencil_dims->x = vol_dims.x - 1;
		}
	}
	if(stencil_dims->y > vol_dims.y) {
			if(vol_dims.y % 2) {
			stencil_half_dims->y = vol_dims.y / 2;
			stencil_dims->y = vol_dims.y;
		} else {
			stencil_half_dims->y = vol_dims.y / 2 - 1;
			stencil_dims->y = vol_dims.y - 1;
		}
	}
	if(stencil_dims->z > vol_dims.z) {
			if(vol_dims.z % 2) {
			stencil_half_dims->z = vol_dims.z / 2;
			stencil_dims->z = vol_dims.z;
		} else {
			stencil_half_dims->z = vol_dims.z / 2 - 1;
			stencil_dims->z = vol_dims.z - 1;
		}
	}*/
	// get stencil size (optimization - use symmetry)
	const int stencil_size = stencil_dims->x * stencil_dims->y * stencil_dims->z;
	kde_stencil = std::make_shared<std::vector<float>>();
	kde_stencil->reserve(stencil_size);	
	// compute stencil values and store them
	const double KDE_constant = 1.0 / (2.0 * M_PI * sqrt(2.0 * M_PI));
	const dvec3 scaled_spacing{spacing.x/h, spacing.y/h, spacing.z/h};
	for (int iz = -stencil_half_dims->z; iz <= stencil_half_dims->z; ++iz) {
		for (int iy = -stencil_half_dims->y; iy <= stencil_half_dims->y; ++iy) {
			for (int ix = -stencil_half_dims->x; ix <= stencil_half_dims->x; ++ix) {				
                double r_sq =	pow( ix * scaled_spacing.x, 2) +
								pow( iy * scaled_spacing.y, 2) +
								pow( iz * scaled_spacing.z, 2);
                kde_stencil->push_back( KDE_constant * exp(-0.5 * r_sq) );
			}
		}
	}
	if(use_scaling_prop.get()) {
		const double denom = (pow(h,3) * extrema_size);
		int index = 0;
		for (int iz = 0; iz < stencil_dims->z; ++iz) {
			for (int iy = 0; iy < stencil_dims->y; ++iy) {
				for (int ix = 0; ix < stencil_dims->x; ++ix) {				
					(*kde_stencil)[index++] /= denom;
				}
			}
		}
	}
	LogProcessorInfo("Stencil dimensions: ");
	LogProcessorInfo(*stencil_dims);
}

void SequentialKDE::fastKDE() {
    auto buff_pair = mesh_in_.getData()->findBuffer(BufferType::IndexAttrib);
    const ivec3 vol_dims = volume_in_.getData()->getDimensions();
    const int vol_xy = vol_dims.x * vol_dims.y;
	const int vol_size = vol_xy * vol_dims.z;
	// check that buffer contains data
    if (buff_pair.first == nullptr) {
		LogProcessorError("Error: no index buffer.");
        volume_out_.setData(std::make_shared<Volume>(vol_dims, DataFloat32::get()));
        return;
    }
    // get extrema data
	const int num_extrema = buff_pair.first->getRepresentation<BufferRAM>()->getSize();
    const int* extrema_indices = static_cast<const int*>(buff_pair.first->getRepresentation<BufferRAM>()->getData());
	// make stencil for fast KDE
	const double h = bandwidth_prop.get();
	makeKDEStencil(num_extrema, h); 
	const int stencil_xy = stencil_dims->x * stencil_dims->y;
	const std::shared_ptr<std::vector<float>> kde_stencil_loc = std::make_shared<std::vector<float>>(*kde_stencil);
	const std::shared_ptr<ivec3> stencil_dims_loc = std::make_shared<ivec3>(*stencil_dims);
	// init output volume
	auto KDE_vr_precision = std::make_shared<VolumeRAMPrecision<float>>(vol_dims);
    float* KDE_raw_ptr = KDE_vr_precision->getDataTyped();
	// for computing max value
	double max_val = 0;
	// Explicitly disable dynamic teams
	omp_set_dynamic(0);     
	// Use 8 threads for all consecutive parallel regions
	omp_set_num_threads(16); 
	// iterate over each extreme value and apply the stencil
 	for (int i = 0; i < num_extrema; ++i) {
		// get the 3D grid coordinates of the current extreme value
		const ivec3 ev_grid_coord(	extrema_indices[i] % vol_dims.x,				// x
									(extrema_indices[i] / vol_dims.x) % vol_dims.y,	// y
									extrema_indices[i] / vol_xy	);					// z
		// compute volume and stencil boundaries for current extremum
		// boundaries in z-dimension
		int vol_z_min = ev_grid_coord.z - stencil_half_dims->z;
		int sten_z_min = 0;
		if(vol_z_min < 0) {
			sten_z_min -= vol_z_min;
			vol_z_min = 0;
		}
		int vol_z_max = ev_grid_coord.z + stencil_half_dims->z + 1;
		int sten_z_max = stencil_dims->z;
		if(vol_z_max > vol_dims.z) {
			sten_z_max -= vol_z_max - vol_dims.z;
			vol_z_max = vol_dims.z;
		}
		// boundaries in y-dimension
		int vol_y_min = ev_grid_coord.y - stencil_half_dims->y;
		int sten_y_min = 0;
		if(vol_y_min < 0) {
			sten_y_min -= vol_y_min;
			vol_y_min = 0;
		}
		int vol_y_max = ev_grid_coord.y + stencil_half_dims->y + 1;
		int sten_y_max = stencil_dims->y;
		if(vol_y_max > vol_dims.y) {
			sten_y_max -= vol_y_max - vol_dims.y;
			vol_y_max = vol_dims.y;
		}
		// boundaries limits in x-dimension
		int vol_x_min = ev_grid_coord.x - stencil_half_dims->x;
		int sten_x_min = 0;
		if(vol_x_min < 0) {
			sten_x_min -= vol_x_min;
			vol_x_min = 0;
		}
		int vol_x_max = ev_grid_coord.x + stencil_half_dims->x + 1;
		int sten_x_max = stencil_dims->x;
		if(vol_x_max > vol_dims.x) {
			sten_x_max -= vol_x_max - vol_dims.x;
			vol_x_max = vol_dims.x;
		}
		// apply stencil centered at ev_grid_coord
		#pragma omp parallel shared(vol_dims, max_val, ev_grid_coord, KDE_raw_ptr, kde_stencil_loc, stencil_dims_loc, stencil_xy, vol_xy, vol_z_min, sten_z_min, sten_z_max, vol_y_min, sten_y_min, sten_y_max, vol_x_min, sten_x_min, sten_x_max)
		{
			double local_max = max_val;
			#pragma omp for nowait
			for(int sten_iz = sten_z_min; sten_iz < sten_z_max; ++sten_iz) {	
				int st_z_offset = sten_iz * stencil_xy;
				int z_offset = (vol_z_min + sten_iz - sten_z_min) * vol_xy;
				for(int sten_iy = sten_y_min; sten_iy < sten_y_max; ++sten_iy) {
					int st_zy_offset = st_z_offset + sten_iy * stencil_dims_loc->x;
					int zy_offset = z_offset + (vol_y_min + sten_iy - sten_y_min) * vol_dims.x;
					for(int sten_ix = sten_x_min; sten_ix < sten_x_max; ++sten_ix) {
						int x_offset = vol_x_min + sten_ix - sten_x_min;
						int offset = zy_offset + x_offset;
						KDE_raw_ptr[offset] += kde_stencil_loc->at(st_zy_offset + sten_ix);
						// store local maximum value seen so far
						if(max_val < KDE_raw_ptr[offset]) max_val = KDE_raw_ptr[offset];
					}
				}
			}
			#pragma omp critical
			{
				if(local_max > max_val) {
					max_val = local_max;
				}
			}
		}
	}
	// create out volume and set properties
	std::shared_ptr<Volume> KDE_vol = std::make_shared<Volume>(KDE_vr_precision);
    KDE_vol->setBasis(volume_in_.getData()->getBasis());
	KDE_vol->setOffset(volume_in_.getData()->getOffset());
	KDE_vol->copyMetaDataFrom(*volume_in_.getData());	
	KDE_vol->dataMap_.valueRange = KDE_vol->dataMap_.dataRange = vec2(0, max_val);
	KDE_vol->setModelMatrix(volume_in_.getData()->getModelMatrix());
	KDE_vol->setWorldMatrix(volume_in_.getData()->getWorldMatrix());
	// set volume on out port
    volume_out_.setData(KDE_vol);
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
