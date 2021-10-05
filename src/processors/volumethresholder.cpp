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

#include <KTH/kerneldensityestimation/processors/volumethresholder.h>

#include <inviwo/core/datastructures/volume/volume.h>
#include <inviwo/core/datastructures/volume/volumeramprecision.h>

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo VolumeThresholder::processorInfo_{
    "org.inviwo.VolumeThresholder",      // Class identifier
    "Volume Thresholder",                // Display name
    "Undefined",              // Category
    CodeState::Experimental,  // Code state
    Tags::None,               // Tags
};
const ProcessorInfo VolumeThresholder::getProcessorInfo() const { return processorInfo_; }

VolumeThresholder::VolumeThresholder()
    : Processor()
    , volume_in_("inport")
	, volume_out_("outport")
    , threshold_prop("threshold","threshold", 0.0001, 0.0, 3, 0.00001) {

    addPort(volume_in_);
	addPort(volume_out_);
    addProperty(threshold_prop);
}

void VolumeThresholder::process() {
    const size3_t vol_dims = volume_in_.getData()->getDimensions();
	const size_t vol_size = vol_dims.x * vol_dims.y * vol_dims.z;
	//const float* in_vol_data = (float*) volume_in_.getData()->getRepresentation<VolumeRAM>();
    auto out_vol = std::make_shared<Volume>(vol_dims, DataFloat32::get());
    const double threshold = threshold_prop.get();	
	double max_val = 0.0;
    out_vol->getEditableRepresentation<VolumeRAM>()->dispatch<void, dispatching::filter::Scalars>([&](auto binary_vol) {
        using BinaryType = util::PrecisionValueType<decltype(binary_vol)>;
        BinaryType* binary_data = binary_vol->getDataTyped();
		const BinaryType* in_vol_data = (const BinaryType*) volume_in_.getData()->getRepresentation<VolumeRAM>()->getData();
		int index = 0;
        for (int iz = 0; iz < vol_dims.z; ++iz) {
            for (int iy = 0; iy < vol_dims.y; ++iy) {
                for (int ix = 0; ix < vol_dims.x; ++ix) {
                    binary_data[index] = in_vol_data[index] < threshold ? 0.0 : 1.0;
					if(max_val < in_vol_data[index]) max_val = in_vol_data[index];
					++index;
				}
            }
        }
    });
	threshold_prop.setMaxValue(max_val);
	threshold_prop.setIncrement(max_val/1000.0);
	out_vol->copyMetaDataFrom(*volume_in_.getData());
	out_vol->dataMap_ = volume_in_.getData()->dataMap_;
    out_vol->dataMap_.dataRange = vec2(0,1);
	out_vol->dataMap_.valueRange = vec2(0,1);
	out_vol->setBasis(volume_in_.getData()->getBasis());
    out_vol->setOffset(volume_in_.getData()->getOffset());
    out_vol->setWorldMatrix(volume_in_.getData()->getWorldMatrix());
    out_vol->setModelMatrix(volume_in_.getData()->getModelMatrix());
	volume_out_.setData(out_vol);
}

}  // namespace inviwo
