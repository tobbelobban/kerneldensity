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

#include <KTH/kerneldensityestimation/processors/volumechannelselector.h>

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo VolumeChannelSelector::processorInfo_{
    "org.inviwo.VolumeChannelSelector",      // Class identifier
    "Volume Channel Selector",                // Display name
    "Undefined",              // Category
    CodeState::Experimental,  // Code state
    Tags::None,               // Tags
};
const ProcessorInfo VolumeChannelSelector::getProcessorInfo() const { return processorInfo_; }

VolumeChannelSelector::VolumeChannelSelector()
    : Processor()
    , volume_in_("Volume_inport")
	, volume_out_("Volume_outport")
	, dimension_("dimension", "Dimension",
              {{"x_dim", "X", Dimension::X},
               {"y_dim", "Y", Dimension::Y},
               {"z_dim", "Z", Dimension::Z}}) {

    addPort(volume_in_);
	addPort(volume_out_);

	addProperty(dimension_);

}

void VolumeChannelSelector::process() {
    const std::shared_ptr<const Volume> volume_in_ptr = volume_in_.getData();
	float min, max;
	std::shared_ptr<Volume> single_channel_vol = volume_in_ptr->getRepresentation<VolumeRAM>()->dispatch<std::shared_ptr<Volume>, dispatching::filter::Float3s>(
		[&](auto vol_in_ram) {
			// set vol dims and size
			const size3_t dims = vol_in_ram->getDimensions();	
			const size_t vol_size = dims.x * dims.y * dims.z;
			const size_t c = static_cast<int>(dimension_.get());
			// get raw ptr to input volume
			const vec3* raw_in_ptr = static_cast<const vec3*>(vol_in_ram->getData());
			// init out vol
			auto out_vol_RAM_repr = std::make_shared<VolumeRAMPrecision<float>>(dims);
			float* out_vol_raw_ptr = out_vol_RAM_repr->getDataTyped();
			// for setting value and data ranges 
			min = std::numeric_limits<float>::max();
			max = std::numeric_limits<float>::min();
			for(int i = 0; i < vol_size; ++i) {
				float res =  raw_in_ptr[i][c];
				out_vol_raw_ptr[i] = res;
				if(res < min) min = res;
				if(res > max) max = res;
			}			
			return std::make_shared<Volume>(out_vol_RAM_repr);
		}	
	);
	
    single_channel_vol->setBasis(volume_in_ptr->getBasis());
	single_channel_vol->setOffset(volume_in_ptr->getOffset());
	single_channel_vol->copyMetaDataFrom(*volume_in_ptr);	
	single_channel_vol->dataMap_.valueRange = single_channel_vol->dataMap_.dataRange = vec2(min, max);
	single_channel_vol->setModelMatrix(volume_in_ptr->getModelMatrix());
	single_channel_vol->setWorldMatrix(volume_in_ptr->getWorldMatrix());
    volume_out_.setData(single_channel_vol);
}

}  // namespace inviwo
