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

#include <KTH/kerneldensityestimation/processors/volumeextrema.h>

#include <inviwo/core/util/glm.h>
#include <inviwo/core/datastructures/volume/volume.h>
#include <inviwo/core/datastructures/volume/volumeramprecision.h>

#include <inviwo/core/datastructures/buffer/buffer.h>
#include <inviwo/core/datastructures/buffer/bufferramprecision.h>
#include <inviwo/core/datastructures/buffer/bufferram.h>

#include <inviwo/core/util/indexmapper.h>
#include <inviwo/core/datastructures/coordinatetransformer.h>

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo VolumeExtrema::processorInfo_{
    "org.inviwo.VolumeExtrema",      // Class identifier
    "Volume Extrema",                // Display name
    "KDE",              // Category
    CodeState::Experimental,  // Code state
    Tags::None,               // Tags
};
const ProcessorInfo VolumeExtrema::getProcessorInfo() const { return processorInfo_; }

VolumeExtrema::VolumeExtrema()
    : Processor()
    , volume_in_("volumeInport")
	, mesh_out_("meshOutport")
	, select_maxima_prop("maxima_selector", "Maxima", true)
    , select_minima_prop("minima_selector", "Minima", false)
    , select_extrema_compare26("tighter_extrema", "Tight bound", false)
	, select_use_abs("use_absolute_comp", "Absolute", false)
	, extrema_count_prop("extrema_count", "Extrema Count", 0, 0, 0)
	{

    addPort(volume_in_);
    addPort(mesh_out_);

	addProperty(select_maxima_prop);
    addProperty(select_minima_prop);
	addProperty(select_extrema_compare26);
	addProperty(select_use_abs);
	addProperty(extrema_count_prop);
}

void VolumeExtrema::process() {
	auto start = std::chrono::high_resolution_clock::now();
	
	const size3_t vol_dims = volume_in_.getData()->getDimensions();
    const size_t vol_size = vol_dims.x * vol_dims.y * vol_dims.z;
    
	const auto world_matrix = volume_in_.getData()->getCoordinateTransformer().getIndexToWorldMatrix();
	const vec3 offset = volume_in_.getData()->getOffset();	
	const auto basis = volume_in_.getData()->getBasis();

	const int want_minima = select_minima_prop.get();
    const int want_maxima = select_maxima_prop.get();
    const bool want_all = want_minima && want_maxima;
	const bool use_abs = select_use_abs.get();
	const bool use_tight_bound = select_extrema_compare26.get();

	const int init_buff_sz = vol_size * 0.005f; // 0.5 % of total volume size

	std::vector<vec3> pos_v;
    std::vector<int> index_v;
	pos_v.reserve(init_buff_sz);
	index_v.reserve(init_buff_sz);

    volume_in_.getData()->getRepresentation<VolumeRAM>()->dispatch<void, dispatching::filter::Scalars>(
		[&](auto vol_pr) {
			
			using VolType = util::PrecisionValueType<decltype(vol_pr)>;
			const VolType* vol_data = vol_pr->getDataTyped();

			// tight bound - compare with 26 adjacent points
			auto check_extrema_26 = [&](const int index, const vec3& coords, const VolType* vol_data)->int {
				const int xy_ = vol_dims.x * vol_dims.y;
				int maxb = 1, minb = 1;
				const VolType curr_val = vol_data[index];
				int z_offset, yz_offset;
				for(int iz = -1; iz < 2; ++iz) {
					if( (coords.z + iz) < 0 || (coords.z + iz) > vol_dims.z-1) continue;
					z_offset = index + iz * xy_;
					for (int iy = -1; iy < 2; ++iy) {			
						if ( (coords.y + iy) < 0 || (coords.y + iy) > vol_dims.y-1) continue;
						yz_offset = z_offset + iy * vol_dims.x;
						for (int ix = -1; ix < 2; ++ix) {
							if ( (coords.x + ix) < 0 || (coords.x + ix) > vol_dims.x-1) continue;
							if( (ix == 0) && (iy == 0) && (iz == 0) ) continue;
							//if( (ix == 0) && (iz == 0) ) continue;
							if(use_abs) {
								maxb = maxb && ( abs((float)curr_val) > abs((float)vol_data[yz_offset + ix]) );				
								minb = minb && ( abs((float)curr_val) < abs((float)vol_data[yz_offset + ix]) );	
							} else {
								maxb = maxb && (curr_val >= vol_data[yz_offset + ix]);				
								minb = minb && (curr_val <= vol_data[yz_offset + ix]);	
							}
							if(!(minb || maxb)) return 0;
						}
					}
				}
				return minb*-1 + maxb;
			};
			
			// compare with 6 adjacent points on grid axes
			auto check_extrema_6 = [&] (const int index, const vec3& coords, const VolType* vol_data)->int {
				const int xy_ = vol_dims.x * vol_dims.y;
				const VolType curr_val = vol_data[index];
				int minb = 1, maxb = 1;
				if(use_abs) {
					minb = minb && (coords.x > 0 ? abs((float)curr_val) <= abs((float)vol_data[index - 1]) : true);
					maxb = maxb && (coords.x > 0 ? abs((float)curr_val) >= abs((float)vol_data[index - 1]) : true);
				
					minb = minb && (coords.x < vol_dims.x - 1 ? abs((float)curr_val) <= abs((float)vol_data[index + 1]) : true);
					maxb = maxb && (coords.x < vol_dims.x - 1 ? abs((float)curr_val) >= abs((float)vol_data[index + 1]) : true);
				
					minb = minb && (coords.y > 0 ? abs((float)curr_val) <= abs((float)vol_data[index - vol_dims.x]) : true);
					maxb = maxb && (coords.y > 0 ? abs((float)curr_val) >= abs((float)vol_data[index - vol_dims.x]) : true);
				
					minb = minb && (coords.y < vol_dims.y - 1 ? abs((float)curr_val) <= abs((float)vol_data[index + vol_dims.x]) : true);
					maxb = maxb && (coords.y < vol_dims.y - 1 ? abs((float)curr_val) >= abs((float)vol_data[index + vol_dims.x]) : true);

					minb = minb && (coords.z > 0 ? abs((float)curr_val) <= abs((float)vol_data[index - xy_]) : true);
					maxb = maxb && (coords.z > 0 ? abs((float)curr_val) >= abs((float)vol_data[index - xy_]) : true);
				
					minb = minb && (coords.z < vol_dims.z - 1 ? abs((float)curr_val) <= abs((float)vol_data[index + xy_]) : true);
					maxb = maxb && (coords.z < vol_dims.z - 1 ? abs((float)curr_val) >= abs((float)vol_data[index + xy_]) : true);
				} else {
					//minb = minb && (coords.x > 0 ? curr_val <= vol_data[index - 1] : true);
					//maxb = maxb && (coords.x > 0 ? curr_val >= vol_data[index - 1] : true);
				
					//minb = minb && (coords.x < vol_dims.x - 1 ? curr_val <= vol_data[index + 1] : true);
					//maxb = maxb && (coords.x < vol_dims.x - 1 ? curr_val >= vol_data[index + 1] : true);
				
					//minb = minb && (coords.y > 0 ? curr_val <= vol_data[index - vol_dims.x] : true);
					//maxb = maxb && (coords.y > 0 ? curr_val >= vol_data[index - vol_dims.x] : true);
				
					//minb = minb && (coords.y < vol_dims.y - 1 ? curr_val <= vol_data[index + vol_dims.x] : true);
					//maxb = maxb && (coords.y < vol_dims.y - 1 ? curr_val >= vol_data[index + vol_dims.x] : true);

					minb = minb && (coords.z > 0 ? curr_val <= vol_data[index - xy_] : true);
					maxb = maxb && (coords.z > 0 ? curr_val >= vol_data[index - xy_] : true);
				
					minb = minb && (coords.z < vol_dims.z - 1 ? curr_val <= vol_data[index + xy_] : true);
					maxb = maxb && (coords.z < vol_dims.z - 1 ? curr_val >= vol_data[index + xy_] : true);
				}
				return minb*-1 + maxb;
			};
			
			// iterate over input volume and find extrema
			int index = 0, min_count = 0, max_count = 0;
			for (int iz = 0; iz < vol_dims.z; ++iz) {
				for (int iy = 0; iy < vol_dims.y; ++iy) {
					for (int ix = 0; ix < vol_dims.x; ++ix) {	
						auto physical_pos = offset + vec3(ix,iy,iz) * vec3(basis[0][0]/(vol_dims.x-1), basis[1][1]/(vol_dims.y-1), basis[2][2]/(vol_dims.z-1));
						if(use_tight_bound) {
							switch (check_extrema_26(index, {ix,iy,iz}, vol_data)) {
								case 0:		// not extrema
									break;
								case -1:	// minima
									min_count++;
									if(want_minima) {
										pos_v.push_back(physical_pos);
										index_v.push_back(index);
									} 
									break;
								case 1:		// maxima
									max_count++;
									if(want_maxima) {
										pos_v.push_back(physical_pos);
										index_v.push_back(index);
									} 
									break;
							}
						} else {
							switch (check_extrema_6(index, {ix,iy,iz}, vol_data)) {
								case 0:
									break;
								case -1:
									min_count++;
									if(want_minima) {
										pos_v.push_back(physical_pos);
										index_v.push_back(index);
									} 
									break;
								case 1:
									max_count++;
									if(want_maxima) {
										pos_v.push_back(physical_pos);
										index_v.push_back(index);
									} 
									break;
							}
						}
						++index;
					}
				}
			}
			
			nr_minima = min_count;
			nr_maxima = max_count;
		}
	);
	
	extrema_count_prop.setMaxValue((long)(want_minima ? nr_minima : 0) + (want_maxima ? nr_maxima : 0));

    auto pos_buffer = std::make_shared<Buffer<vec3>>(std::make_shared<BufferRAMPrecision<vec3>>(pos_v));
    auto index_buffer = std::make_shared<Buffer<int>>(std::make_shared<BufferRAMPrecision<int>>(index_v));

    auto stop = std::chrono::high_resolution_clock::now();

	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::string time_str = "Create extrema buffer: " + std::to_string(duration.count()) + " ms";
	std::string minima_str = "Nr minima: " + std::to_string(nr_minima);
	std::string maxima_str = "Nr maxima: " + std::to_string(nr_maxima);
    LogProcessorInfo(time_str);
	
	if(want_minima) LogProcessorInfo(minima_str);
	if(want_maxima) LogProcessorInfo(maxima_str);

	auto mesh = std::make_shared<Mesh>(DrawType::Points, ConnectivityType::None);
	mesh->addBuffer(BufferType::PositionAttrib, pos_buffer);
	mesh->addBuffer(BufferType::IndexAttrib, index_buffer);

	mesh_out_.setData(mesh);
}

}  // namespace inviwo
