
####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was ViennaHRLEConfig.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

option(OMP_PARALLELIZE "Use OpenMP to parallelize computations." ON)

###############################################
#  compiler dependent settings for OpenMP
###############################################
if(OMP_PARALLELIZE)
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp=libomp")
  elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
  elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopenmp")
  elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")
  endif()
endif()

# set different include path for local builds
if(NOT ViennaHRLE_COMMON_TARGET)
  include("${CMAKE_CURRENT_LIST_DIR}/ViennaHRLETargets.cmake")
endif()

SET(VIENNAHRLE_INCLUDE_DIRS "/home/xaver/Code/ViennaLS/Simulations/GeometricBosch/ViennaTools/viennahrle/install/ViennaHRLE/include")

check_required_components("ViennaHRLE")
