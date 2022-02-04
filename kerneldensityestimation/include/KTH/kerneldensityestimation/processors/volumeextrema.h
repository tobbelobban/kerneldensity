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

#include <KTH/kerneldensityestimation/kerneldensityestimationmoduledefine.h>
#include <inviwo/core/processors/processor.h>
#include <inviwo/core/properties/boolproperty.h>
#include <inviwo/core/properties/ordinalproperty.h>
#include <inviwo/core/ports/volumeport.h>
#include <inviwo/core/ports/bufferport.h>
#include <inviwo/core/ports/meshport.h>

namespace inviwo {

class IVW_MODULE_KERNELDENSITYESTIMATION_API VolumeExtrema : public Processor {
public:
	VolumeExtrema();
	virtual ~VolumeExtrema() = default;

	virtual void process() override;

	virtual const ProcessorInfo getProcessorInfo() const override;
	static const ProcessorInfo processorInfo_;

private:
	// ports
	VolumeInport volume_in_;
	MeshOutport mesh_out_;

	// properties
	BoolProperty select_maxima;		// if user wants maxima
	BoolProperty select_minima;		// if user wants minima
	BoolProperty select_use_N26;		// if user wants to use N_26 neighbourhood when comparing vertices
	BoolProperty select_use_abs;		// if user wants to use absolute values during comparison

	size_t nr_maxima = 0;
	size_t nr_minima = 0;

	bool use_abs;

	// functions
	int extreme_value_check_N26(	const size_t index, 
					const size3_t coords,
					const size3_t vol_dims,
					const float* vol_data	);
	
	int extreme_value_check_N6(	const size_t index,
					const size3_t coords,
					const size3_t vol_dims,
					const float* vol_data	);

};

	

}  // namespace inviwo
