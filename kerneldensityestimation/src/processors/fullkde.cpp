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

#include <inviwo/kerneldensityestimation/processors/fullkde.h>

#include <inviwo/core/util/glm.h>
#include <inviwo/core/datastructures/volume/volume.h>
#include <inviwo/core/datastructures/volume/volumeramprecision.h>

#include <inviwo/core/datastructures/buffer/buffer.h>
#include <inviwo/core/datastructures/buffer/bufferramprecision.h>
#include <inviwo/core/datastructures/buffer/bufferram.h>

#include <inviwo/core/util/indexmapper.h>
#include <inviwo/core/util/threadpool.h>

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo FullKDE::processorInfo_{
    "org.inviwo.FullKDE",      // Class identifier
    "Full KDE",                // Display name
    "Undefined",              // Category
    CodeState::Experimental,  // Code state
    Tags::None,               // Tags
};
const ProcessorInfo FullKDE::getProcessorInfo() const { return processorInfo_; }

FullKDE::FullKDE()
    : Processor()
    , volInport("inport")
	, volOutport("outport") {

    addPort(volInport);
    addPort(volOutport);
}

void FullKDE::process() {	
	auto vol = 
        volInport.getData()->getRepresentation<VolumeRAM>()->dispatch<std::shared_ptr<Volume>>(
            [&](auto volInPrecision) {
                //using ValType = util::PrecisionValueType<decltype(volInPrecision)>();
                const size3_t dims = volInPrecision->getDimensions();
                auto volInData = (const float*)volInPrecision->getDataTyped();
                int index = 0;
                std::vector<int> extrema; // currently local maxima only
                int extremaCount = 0; 
                int nonExtremaCount = 0;
                for (int iz = 0; iz < dims.z; ++iz) {
                    for (int iy = 0; iy < dims.y; ++iy) {
                        for (int ix = 0; ix < dims.x; ++ix) {
                            bool isMaxima = true, isMinima = true;
                            float currVal = volInData[index];
                            if (ix != 0 && volInData[index - 1] > currVal) isMaxima = false;
                            if (ix != 0 && volInData[index - 1] < currVal) isMinima = false;

                            if (ix != dims.x - 1 && volInData[index + 1] > currVal) isMaxima = false;
                            if (ix != dims.x - 1 && volInData[index + 1] < currVal) isMinima = false;

							if (iy != 0 && volInData[index - dims.x] > currVal) isMaxima = false;
                            if (iy != 0 && volInData[index - dims.x] < currVal) isMinima = false;

                            if (iy != dims.y - 1 && volInData[index + dims.x] > currVal) isMaxima = false;
							if (iy != dims.y - 1 && volInData[index + dims.x] < currVal) isMinima = false;

                            if (iz != 0 && volInData[index - dims.x*dims.y] > currVal) isMaxima = false;
							if (iz != 0 && volInData[index - dims.x*dims.y] < currVal) isMinima = false;

                            if (iz != dims.z - 1 && volInData[index + dims.x*dims.y] > currVal) isMaxima = false;
							if (iz != dims.z - 1 && volInData[index + dims.x*dims.y] < currVal) isMinima = false;

                        if (isMaxima != isMinima) {
                                extrema.push_back(index);
                            extremaCount += 1;
                        } else {
                            nonExtremaCount += 1;
                        }
                        ++index;
                    }
                }
            }
			
			IVW_ASSERT(extremaCount + nonExtremaCount == dims.x * dims.y * dims.z, "Error! Incorrect extrema!");
			// volume stored in x -> y -> z
            /*LogProcessorInfo(volInPrecision->getAsDouble({0,0,1}));
            LogProcessorInfo(volInPrecision->getAsDouble({1,0,1}));
            LogProcessorInfo(volInPrecision->getAsDouble({2,0,1}));
            LogProcessorInfo(volInPrecision->getAsDouble({3,0,1}));

			LogProcessorInfo(volInData[dims.x*dims.y]);
            LogProcessorInfo(volInData[dims.x * dims.y + 1]);
            LogProcessorInfo(volInData[dims.x * dims.y + 2]);
            LogProcessorInfo(volInData[dims.x * dims.y + 3]);*/

            auto volOutPrecision = std::make_shared<VolumeRAMPrecision<float>>(dims);
            float* volOutData = volOutPrecision->getDataTyped();
            double h = 0.009;
            index = 0;                      
            for (int iz = 0; iz < dims.z; ++iz) {
                for (int iy = 0; iy < dims.y; ++iy) {
                    for (int ix = 0; ix < dims.x; ++ix) {
                        dvec3 worldPos = {ix * 0.005, iy * 0.001, iz * 0.005};
						double val = 0;
                        double radius_sq;
						for (int i = 0; i < extrema.size(); ++i) {
							dvec3 extremaWorldPos = {(extrema[i] % dims.x) * 0.005,
													(extrema[i] / dims.x) % dims.y * 0.001,
													(int)(extrema[i] / (dims.x * dims.y)) * 0.005};
							dvec3 scaledPos = (worldPos - extremaWorldPos);
							scaledPos /= h;
							radius_sq = scaledPos.x * scaledPos.x + scaledPos.y * scaledPos.y + scaledPos.z * scaledPos.z;
							val += (1.0 / (2.0 * M_PI * sqrt(2.0 * M_PI))) * exp(-0.5 * radius_sq);
						}
						val /= (h * extrema.size());              
                        volOutData[index] = (float)val;
                        ++index;
                    }
                }
            }
			return std::make_shared<Volume>(volOutPrecision);
        }
	);
	vol->copyMetaDataFrom(*volInport.getData());
    vol->dataMap_ = volInport.getData()->dataMap_;
	//vol->dataMap_.valueRange = vec2(0, maxVal);
    vol->setBasis(volInport.getData()->getBasis());
    vol->setOffset(volInport.getData()->getOffset());
	vol->setWorldMatrix(volInport.getData()->getWorldMatrix());
    vol->setModelMatrix(volInport.getData()->getModelMatrix());
    volOutport.setData(vol);
}

}  // namespace inviwo
