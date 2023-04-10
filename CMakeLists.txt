cmake_minimum_required(VERSION 3.18)
project(_correa VERSION 0.1.0 LANGUAGES CXX)
find_package(Python 3.8 COMPONENTS Interpreter Development.Module REQUIRED)
#set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# C++ 17 required
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_CXX_EXTENSIONS ON)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

option(CORREA_BUILD_PYTHON_BINDINGS OFF)
option(CORREA_BUILD_EXECUTABLES OFF)

# Default to Release
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/ext/nanobind)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/ext/hera)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/ext/hera/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include/Clustering/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include/Curvature/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include/DataStruct/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include/InOut/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include/Metrics/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include/OT1/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include/PH0/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)



if(CORREA_BUILD_PYTHON_BINDINGS)
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/bindings)
  nanobind_add_module(_correa src/correa.cpp)
  target_link_libraries(_correa PRIVATE hera)
endif()


if(CORREA_BUILD_EXECUTABLES)
    add_executable(2DShape src/2DShape.cpp)
    add_executable(Comp2DShape src/Comp2DShapes.cpp)
endif()
