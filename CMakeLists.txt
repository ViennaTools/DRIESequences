cmake_minimum_required(VERSION 3.4)

project("GeometricBosch")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

## Enable Clang sanitizer for debug builds
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address -fsanitize=thread -fsanitize=memory" CACHE STRING "")
  set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUGS} -fno-omit-frame-pointer -fsanitize=address -fsanitize=thread -fsanitize=memory" CACHE STRING "")
endif()

# automatically append standard path for ViennaHRLEConfig.cmake
# might need to be changed for different install directories
# not necessary when installed to /usr/local (not recommended)
SET(ViennaToolsDir "ViennaTools/")
SET(ViennaHRLE_DIR "${ViennaToolsDir}/viennahrle/install/ViennaHRLE/lib/cmake/ViennaHRLE")
SET(ViennaLS_DIR "${ViennaToolsDir}/viennals/install/ViennaLS/lib/cmake/ViennaLS")

find_package(ViennaHRLE REQUIRED)
find_package(ViennaLS REQUIRED)

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