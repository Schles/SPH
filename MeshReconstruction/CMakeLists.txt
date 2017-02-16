# Use static anttweakbar.
add_definitions(-DTW_NO_LIB_PRAGMA -DTW_STATIC)

# Use static GLEW.
add_definitions(-DGLEW_STATIC)

# Use static freeglut.
add_definitions(-DFREEGLUT_LIB_PRAGMAS=0 -DFREEGLUT_STATIC)

# CompactNSearch
add_definitions(-D_USE_MATH_DEFINES)


if (UNIX)
find_package(OpenMP)
if (OPENMP_FOUND)
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif(OPENMP_FOUND)
endif (UNIX)

file(TO_NATIVE_PATH "${PROJECT_SOURCE_DIR}/resources/" GLVIZ_RESOURCES_DIR)
configure_file(config.hpp.in "${CMAKE_CURRENT_BINARY_DIR}/config.hpp")

file(GLOB SOURCES "kernel/*.cpp" "ex/*.cpp" "kernel/*.h" "ex/*.h")

add_executable(FluidMeshReconstruction
    main.cpp
    "../src/kernel/Kernel.h"
    "../src/kernel/Kernel.cpp"
    ${SOURCES}
)

include_directories(
    ${CMAKE_CURRENT_BINARY_DIR}
    ${glviz_INCLUDE_DIRS}
    "../CompactNSearch/include"
    "../src/kernel"
)

add_dependencies(FluidMeshReconstruction
    glviz
    CompactNSearch
)

target_link_libraries(FluidMeshReconstruction
    ${glviz_LIBRARIES}
    CompactNSearch
)
