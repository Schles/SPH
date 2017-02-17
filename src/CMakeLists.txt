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

file(GLOB SOURCES "kernel/*.cpp" "sph/*.cpp" "kernel/*.h"  "sph/*.h" "particle/*.cpp" "particle/*.h")

add_executable(FluidSim
    main.cpp
    Scenes.cpp
    ${SOURCES}
)

include_directories(
    ${CMAKE_CURRENT_BINARY_DIR}
    ${glviz_INCLUDE_DIRS}
    "../CompactNSearch/include"
    "kernel"
    "util"
    "sph"	
    "particle"
)

add_dependencies(FluidSim
    glviz
    CompactNSearch
)

target_link_libraries(FluidSim
    ${glviz_LIBRARIES}
    CompactNSearch
)
