cmake_minimum_required(VERSION 3.18)
project(correa VERSION 1.0.0 LANGUAGES CXX)
find_package(Python 3.8 COMPONENTS Interpreter Development.Module REQUIRED)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# C++ 17 required
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_CXX_EXTENSIONS ON)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
option(CORREA_BUILD_EXECUTABLES "build executables" OFF)	
option(CORREA_BUILD_PYTHON_BINDINGS "build python bindings" ON)
# Boost

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
find_package(BLAS)
find_package(LAPACK)


set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include/Clustering/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include/Curvature/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include/DataStruct/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include/InOut/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include/Metrics/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include/OT1/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include/PH0/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/ext/hera/include)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/ext/hera)


add_library(correa INTERFACE)
if (CORREA_BUILD_EXECUTABLES)
	add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/src)
endif()

if(CORREA_BUILD_PYTHON_BINDINGS)
	#add_subdirectory(bindings)
	add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/ext/nanobind)
	include_directories(${CMAKE_CURRENT_SOURCE_DIR}/bindings)
	nanobind_add_module(_correa src/correa.cpp)
	target_compile_definitions(_correa PRIVATE CORREA_PYTHON_BINDINGS)
	target_link_libraries(_correa PRIVATE ${libraries} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} hera)
	file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/bindings/correa/__init__.py DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/build/correa)
	# Copy the compiled module to build/correa after it's built
	add_custom_command(TARGET _correa POST_BUILD
		COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:_correa> ${CMAKE_CURRENT_SOURCE_DIR}/build/correa/
	)
endif()
