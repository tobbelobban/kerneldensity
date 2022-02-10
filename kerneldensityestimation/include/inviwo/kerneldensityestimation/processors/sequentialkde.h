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

#pragma once

#include <inviwo/kerneldensityestimation/kerneldensityestimationmoduledefine.h>
#include <inviwo/core/processors/processor.h>
#include <inviwo/core/properties/boolproperty.h>
#include <inviwo/core/properties/ordinalproperty.h>
#include <inviwo/core/ports/volumeport.h>
#include <inviwo/core/ports/meshport.h>

namespace inviwo {

class IVW_MODULE_KERNELDENSITYESTIMATION_API SequentialKDE : public Processor {
public:
	SequentialKDE();
	virtual ~SequentialKDE() = default;

	virtual void process() override;

	virtual const ProcessorInfo getProcessorInfo() const override;
	static const ProcessorInfo processorInfo_;

private:
    	// ports
	VolumeInport volume_in_;
	MeshInport mesh_in_;
	VolumeOutport volume_out_;

	// properties
	FloatProperty bandwidth_prop;		// for choosing KDE bandwidth
	FloatProperty cutoff_prop;		// for choosing stencil size in fastKDE()								                              
	BoolProperty fast_KDE_prop;		// if user wants fastKDE
	BoolProperty use_scaling_prop;		// if user wants to scale result

	// containers for stencil
	std::shared_ptr<std::vector<float>> kde_stencil;
	std::shared_ptr<ivec3> stencil_dims;
	std::shared_ptr<ivec3> stencil_half_dims;

	// functions 
	void KDE();
	void fastKDE();
	void makeKDEStencil(const size_t nr_extrema, const double h);
	
};

}  // namespace inviwo
