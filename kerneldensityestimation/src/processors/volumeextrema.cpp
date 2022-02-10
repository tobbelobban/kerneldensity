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

#include <inviwo/kerneldensityestimation/processors/volumeextrema.h>

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
    , select_maxima("maxima_selector", "Maxima", true)
    , select_minima("minima_selector", "Minima", false)
    , select_use_N26("n26", "Use N_26 neighbourhood", false)
{

    addPort(volume_in_);
    addPort(mesh_out_);

    addProperty(select_maxima);
    addProperty(select_minima);
    addProperty(select_use_N26);
}

/*
Check if vertex at coords is an extrema in vol_data by comparing with
neighbourhood consisting of vertices that share either an edge, face
or voxel with vertex at coords. In total 26 vertices.
Returns:
0 - not extrema
-1 - minima
1 - maxima
*/
int VolumeExtrema::extreme_value_check_N26(const size_t index,
    const size3_t coords,
    const size3_t vol_dims,
    const float* vol_data)
{
    int maxb = 1, minb = 1;
    const size_t d_x = vol_dims.x;
    const size_t d_xy = d_x * vol_dims.y;
    const size_t vol_size = vol_dims.x * vol_dims.y * vol_dims.z;
    const float curr_val = vol_data[index];
    size_t z_offset, yz_offset;
    for (int iz = -1; iz < 2; ++iz) {
        if ((coords.z + iz) < 0 || (coords.z + iz) >= vol_dims.z) continue;
        z_offset = index + iz * d_xy;
        for (int iy = -1; iy < 2; ++iy) {
            if ((coords.y + iy) < 0 || (coords.y + iy) >= vol_dims.y) continue;
            yz_offset = z_offset + iy * d_x;
            for (int ix = -1; ix < 2; ++ix) {
                if ((coords.x + ix) < 0 || (coords.x + ix) >= vol_dims.x) continue;
                if ((ix == 0) && (iy == 0) && (iz == 0)) continue;
                maxb = maxb && (curr_val > vol_data[yz_offset + ix]);
                minb = minb && (curr_val < vol_data[yz_offset + ix]);
                if (!(minb || maxb)) return 0;
            }
        }
    }
    return minb * -1 + maxb;
}

/*
Check if vertex at coords is an extrema in vol_data by comparing with
neighbourhood consisting of vertices that share an edge with vertex
at coords. In total 6 vertices.
Returns:
0 - not extrema
-1 - minima
1 - maxima
*/
int VolumeExtrema::extreme_value_check_N6(const size_t index,
    const size3_t coords,
    const size3_t vol_dims,
    const float* vol_data)
{
    int maxb = 1, minb = 1;
    const size_t d_x = vol_dims.x;
    const size_t d_xy = d_x * vol_dims.y;
    const size_t vol_size = vol_dims.x * vol_dims.y * vol_dims.z;
    // out of bounds
    if (index > vol_size) return 0;
    const float curr_val = vol_data[index];
    // compare in x dim
    for (int ix = -1; ix < 2; ix += 2) {
        if ((coords.x + ix) < 0 || (coords.x + ix) >= vol_dims.x) continue;
        maxb = maxb && (curr_val > vol_data[index + ix]);
        minb = minb && (curr_val < vol_data[index + ix]);
    }
    // compare in y dim
    for (int iy = -1; iy < 2; iy += 2) {
        if ((coords.y + iy) < 0 || (coords.y + iy) >= vol_dims.y) continue;
        maxb = maxb && (curr_val > vol_data[index + iy * d_x]);
        minb = minb && (curr_val < vol_data[index + iy * d_x]);
    }
    // compare in z dim
    for (int iz = -1; iz < 2; iz += 2) {
        if ((coords.z + iz) < 0 || (coords.z + iz) >= vol_dims.z) continue;
        maxb = maxb && (curr_val > vol_data[index + iz * d_xy]);
        minb = minb && (curr_val < vol_data[index + iz * d_xy]);
    }
    return minb * -1 + maxb;
}

void VolumeExtrema::process() {
    // input volume pointer
    // NOTE! assumption: input volume contains single-precision floats
    const std::shared_ptr<const Volume> in_v_ptr = volume_in_.getData();
    if (in_v_ptr->getDataFormat()->getId() != DataFormatId::Float32) return;

    // get access to raw input data
    const VolumeRAM* v_ram = in_v_ptr->getRepresentation<VolumeRAM>();
    const float* in_v_raw_ptr = static_cast<const float*>(v_ram->getData());

    // get input volume properties
    const auto world_matrix = in_v_ptr->getCoordinateTransformer().getIndexToWorldMatrix();
    const vec3 offset = in_v_ptr->getOffset();
    const auto basis = in_v_ptr->getBasis();
    const size3_t vol_dims = in_v_ptr->getDimensions();
    const size_t vol_size = vol_dims.x * vol_dims.y * vol_dims.z;

    // get user options
    const bool want_minima = select_minima.get();
    const bool want_maxima = select_maxima.get();
    const bool use_N26 = select_use_N26.get();

    // init buffers for extrema
    const int init_buff_sz = (int)(vol_size * 0.005f); // 0.5 % of total volume size
    std::vector<vec3> pos_v;
    std::vector<int> index_v;
    pos_v.reserve(init_buff_sz);
    index_v.reserve(init_buff_sz);

    // iterate over input volume and find extrema
    size_t index = 0, num_extrema = 0;
    for (int iz = 0; iz < vol_dims.z; ++iz) {
        for (int iy = 0; iy < vol_dims.y; ++iy) {
            for (int ix = 0; ix < vol_dims.x; ++ix) {

                int res = 0;
                if (use_N26) {
                    res = extreme_value_check_N26(index, { ix,iy,iz }, vol_dims, in_v_raw_ptr);
                }
                else {
                    res = extreme_value_check_N6(index, { ix,iy,iz }, vol_dims, in_v_raw_ptr);
                }

                if (res == 0) continue; //no extremum here

                // get current vertex real world position
                vec3 physical_pos = offset + vec3(ix, iy, iz) * vec3(basis[0][0] / (vol_dims.x - 1),
                    basis[1][1] / (vol_dims.y - 1),
                    basis[2][2] / (vol_dims.z - 1));

                switch (res) {
                case -1:	// minimum						
                    if (want_minima) {
                        pos_v.push_back(physical_pos);
                        index_v.push_back((int)index);
                        ++num_extrema;
                    }
                    break;
                case 1:		// maximum
                    if (want_maxima) {
                        pos_v.push_back(physical_pos);
                        index_v.push_back((int)index);
                        ++num_extrema;
                    }
                    break;
                }
                ++index;
            }
        }
    }

    // print number of extrema
    LogProcessorInfo("Number of extrema: " << num_extrema);

    // create buffers for mesh
    auto pos_buffer = std::make_shared<Buffer<vec3>>(std::make_shared<BufferRAMPrecision<vec3>>(pos_v));
    auto index_buffer = std::make_shared<Buffer<int>>(std::make_shared<BufferRAMPrecision<int>>(index_v));
    auto mesh = std::make_shared<Mesh>(DrawType::Points, ConnectivityType::None);
    mesh->addBuffer(BufferType::PositionAttrib, pos_buffer);
    mesh->addBuffer(BufferType::IndexAttrib, index_buffer);

    // set mesh ptr on outport
    mesh_out_.setData(mesh);
}

}  // namespace inviwo
