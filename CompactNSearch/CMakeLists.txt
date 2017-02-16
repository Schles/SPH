cmake_minimum_required(VERSION 3.2)

project(CompactNSearch)

# Visual studio solution directories.
set_property(GLOBAL PROPERTY USE_FOLDERS on)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (UNIX)
find_package(OpenMP)
if (OPENMP_FOUND)
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif(OPENMP_FOUND)
endif (UNIX)

SET(CMAKE_DEBUG_POSTFIX "_d")

set (HEADER_FILES 
		include/Config.h
		include/CompactNSearch.h
		include/PointSet.h
		include/DataStructures.h)

include_directories("include"
		../extern/tbb/include)

add_library(CompactNSearch
		${HEADER_FILES}
		src/CompactNSearch.cpp
)

if(UNIX)
	add_dependencies(CompactNSearch    tbb )
	target_link_libraries(CompactNSearch tbb)
endif(UNIX)

install(FILES "include/CompactNSearch" ${HEADER_FILES}
	DESTINATION include/)

install(TARGETS CompactNSearch
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib)

option(BUILD_DEMO "Build example of how to use this library."
		ON)
if(BUILD_DEMO)
	add_subdirectory(demo)
endif(BUILD_DEMO)

