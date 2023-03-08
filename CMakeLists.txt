cmake_minimum_required(VERSION 3.18)
project(correa VERSION 0.1.0 LANGUAGES CXX)
find_package(Python 3.8 COMPONENTS Interpreter Development.Module REQUIRED)
#set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS ON)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

# Default to Release

if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/ext/nanobind)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/ext/hera)

include_directories(/opt/homebrew/Cellar/boost/1.81.0_1/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/ext/hera/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/Clustering/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/Curvature/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/DataStruct/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/InOut/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/Metrics/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/OT1/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/PH0/include)

#target_link_libraries(correa PRIVATE hera)

nanobind_add_module(_correa src/correa.cpp)