cmake_minimum_required(VERSION 3.4)

project("GeometricBosch")

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

## Enable Clang sanitizer for debug builds
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address -fsanitize=thread -fsanitize=memory" CACHE STRING "")
  set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUGS} -fno-omit-frame-pointer -fsanitize=address -fsanitize=thread -fsanitize=memory" CACHE STRING "")
endif()

# build dependencies
set(DEPENDENCIES_DIR ${CMAKE_SOURCE_DIR}/dependencies)

# Include all external dependencies
include(ExternalProject)
add_custom_target(buildDependencies)
add_subdirectory(external/upstream)

find_package(ViennaLS  PATHS ${ViennaLS_DIR} QUIET)

if(ViennaLS_FOUND)
  option(BUILD_TEST "Build unit test." OFF)
  if(BUILD_TEST)
    add_executable(Test Test/Test.cpp)
    target_include_directories(Test PUBLIC ${VIENNALS_INCLUDE_DIRS})
    target_link_libraries(Test PRIVATE ${VIENNALS_LIBRARIES})
  endif()

  SET(DEM3D "DEM3D")
  add_executable(${DEM3D} ${DEM3D}.cpp)
  target_include_directories(${DEM3D} PUBLIC ${VIENNALS_INCLUDE_DIRS})
  target_link_libraries(${DEM3D} PRIVATE ${VIENNALS_LIBRARIES})

  set(DEM2D "DEM2D")
  add_executable(${DEM2D} ${DEM2D}.cpp)
  target_include_directories(${DEM2D} PUBLIC ${VIENNALS_INCLUDE_DIRS})
  target_link_libraries(${DEM2D} PRIVATE ${VIENNALS_LIBRARIES})

  SET(DREAM "DREAM")
  add_executable(${DREAM} ${DREAM}.cpp)
  target_include_directories(${DREAM} PUBLIC ${VIENNALS_INCLUDE_DIRS})
  target_link_libraries(${DREAM} PRIVATE ${VIENNALS_LIBRARIES})

  SET(DREM3D "DREM3D")
  add_executable(${DREM3D} ${DREM3D}.cpp)
  target_include_directories(${DREM3D} PUBLIC ${VIENNALS_INCLUDE_DIRS})
  target_link_libraries(${DREM3D} PRIVATE ${VIENNALS_LIBRARIES})
else()
  message(STATUS "ViennaLS could not be found. Please install ViennaLS first by runnig `make buildDependencies` 
  or provide the path with `-D ViennaLS_DIR=/path/to/your/ViennaLS/install`.")      
endif()