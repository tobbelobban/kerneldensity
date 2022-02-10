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

#include <inviwo/kerneldensityestimation/processors/volumesubsetdrawer.h>

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo VolumeSubsetDrawer::processorInfo_{
    "org.inviwo.VolumeSubsetDrawer",      // Class identifier
    "Volume Subset Drawer",                // Display name
    "Undefined",              // Category
    CodeState::Experimental,  // Code state
    Tags::None,               // Tags
};
const ProcessorInfo VolumeSubsetDrawer::getProcessorInfo() const { return processorInfo_; }

VolumeSubsetDrawer::VolumeSubsetDrawer()
    : Processor()
    , volume_inport_("VolumeInport")
	, mesh_outport_("MeshOutport") {

    addPort(volume_inport_);
    addPort(mesh_outport_);
}

void VolumeSubsetDrawer::process() {
	auto vol = volume_inport_.getData();
	auto basis = vol->getBasis();
	auto offset = vol->getOffset();
	auto dims = vol->getDimensions();
	vec3 spacing = vec3(basis[0][0]/(dims.x-1), basis[1][1]/(dims.y-1), basis[2][2]/(dims.z-1));
	using MyMesh = TypedMesh<buffertraits::PositionsBuffer,buffertraits::ColorsBuffer>;
	std::vector<MyMesh::Vertex> vertices;

	vertices.emplace_back(spacing * vec3(590,430,0)	+ offset, vec4(0,0,0,1));
	vertices.emplace_back(spacing * vec3(590,430,500) + offset, vec4(0,0,0,1));

	vertices.emplace_back(spacing * vec3(657,274,0)	+ offset, vec4(0,0,0,1));
	vertices.emplace_back(spacing * vec3(657,274,500)	+ offset, vec4(0,0,0,1));

	vertices.emplace_back(spacing * vec3(596,435,0)	+ offset, vec4(0,0,0,1));
	vertices.emplace_back(spacing * vec3(596,435,500)	+ offset, vec4(0,0,0,1));

	vertices.emplace_back(spacing * vec3(663,279,0)	+ offset, vec4(0,0,0,1));
	vertices.emplace_back(spacing * vec3(663,279,500)	+ offset, vec4(0,0,0,1));

	MyMesh mesh;
	mesh.addVertices(vertices);
	auto ib = mesh.addIndexBuffer(DrawType::Lines, ConnectivityType::None);
	ib->add({1, 0, 2, 3, 7, 5, 4, 6, 2, 0, 4, 5, 1, 3, 7, 6, 0, 4, 1, 5, 2, 6, 3, 7});
	mesh_outport_.setData(std::make_shared<Mesh>(mesh));
}

}  // namespace inviwo
