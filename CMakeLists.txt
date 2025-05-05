cmake_minimum_required(VERSION 3.15)

project("DRIESequences")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set default build type
set(DEFAULT_BUILD_TYPE "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${DEFAULT_BUILD_TYPE}' as none was specified.")
  set(CMAKE_BUILD_TYPE
      "${DEFAULT_BUILD_TYPE}"
      CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel"
                                               "RelWithDebInfo")
endif()

include(cmake/cpm.cmake)

CPMFindPackage(
  NAME ViennaLS
  VERSION 4.3.1
  GIT_REPOSITORY "https://github.com/ViennaTools/ViennaLS")

SET(DEM3D "DEM3D")
add_executable(${DEM3D} ${DEM3D}.cpp)
target_include_directories(${DEM3D} PUBLIC ${VIENNALS_INCLUDE_DIRS})
target_link_libraries(${DEM3D} PRIVATE ViennaTools::ViennaLS)

set(DEM2D "DEM2D")
add_executable(${DEM2D} ${DEM2D}.cpp)
target_include_directories(${DEM2D} PUBLIC ${VIENNALS_INCLUDE_DIRS})
target_link_libraries(${DEM2D} PRIVATE ViennaTools::ViennaLS)

SET(DREAM "DREAM")
add_executable(${DREAM} ${DREAM}.cpp)
target_include_directories(${DREAM} PUBLIC ${VIENNALS_INCLUDE_DIRS})
target_link_libraries(${DREAM} PRIVATE ViennaTools::ViennaLS)

SET(DREM3D "DREM3D")
add_executable(${DREM3D} ${DREM3D}.cpp)
target_include_directories(${DREM3D} PUBLIC ${VIENNALS_INCLUDE_DIRS})
target_link_libraries(${DREM3D} PRIVATE ViennaTools::ViennaLS)
