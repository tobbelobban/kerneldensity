ivw_module(KernelDensityEstimation)

set(HEADER_FILES
    include/KTH/kerneldensityestimation/kerneldensityestimationmodule.h
    include/KTH/kerneldensityestimation/kerneldensityestimationmoduledefine.h
    include/KTH/kerneldensityestimation/processors/fullkde.h
    include/KTH/kerneldensityestimation/processors/fullkdecl.h
    include/KTH/kerneldensityestimation/processors/sequentialkde.h
    include/KTH/kerneldensityestimation/processors/volumeextrema.h
    include/KTH/kerneldensityestimation/processors/volumeextremapool.h
    include/KTH/kerneldensityestimation/processors/volumethresholder.h
)
ivw_group("Header Files" ${HEADER_FILES})

set(SOURCE_FILES
    src/kerneldensityestimationmodule.cpp
    src/processors/fullkde.cpp
    src/processors/fullkdecl.cpp
    src/processors/sequentialkde.cpp
    src/processors/volumeextrema.cpp
    src/processors/volumeextremapool.cpp
    src/processors/volumethresholder.cpp
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
