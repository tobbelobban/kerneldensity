ivw_module(KernelDensityEstimation)

set(HEADER_FILES
    include/inviwo/kerneldensityestimation/kerneldensityestimationmodule.h
    include/inviwo/kerneldensityestimation/kerneldensityestimationmoduledefine.h
    include/inviwo/kerneldensityestimation/processors/fullkde.h
    include/inviwo/kerneldensityestimation/processors/sequentialkde.h
    include/inviwo/kerneldensityestimation/processors/slicesequentialkde.h
    include/inviwo/kerneldensityestimation/processors/volumechannelselector.h
    include/inviwo/kerneldensityestimation/processors/volumeextrema.h
    include/inviwo/kerneldensityestimation/processors/volumesubsetdrawer.h
)
ivw_group("Header Files" ${HEADER_FILES})

set(SOURCE_FILES
    src/kerneldensityestimationmodule.cpp
    src/processors/fullkde.cpp
    src/processors/sequentialkde.cpp
    src/processors/slicesequentialkde.cpp
    src/processors/volumechannelselector.cpp
    src/processors/volumeextrema.cpp
    src/processors/volumesubsetdrawer.cpp
)
ivw_group("Source Files" ${SOURCE_FILES})

set(SHADER_FILES
	# Add shaders
)
ivw_group("Shader Files" ${SHADER_FILES})

set(TEST_FILES
    tests/unittests/kerneldensityestimation-unittest-main.cpp
)
ivw_add_unittest(${TEST_FILES})

ivw_create_module(${SOURCE_FILES} ${HEADER_FILES} ${SHADER_FILES})

# Add shader directory to install package
#ivw_add_to_module_pack(${CMAKE_CURRENT_SOURCE_DIR}/glsl)
