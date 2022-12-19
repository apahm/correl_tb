#[=======================================================================[.rst:
FindIPP
-------

Finds the intel IPP libraries.

Imported Targets
^^^^^^^^^^^^^^^^

This module provides the following imported targets, if found:

``ipp::core``
  The ippcore library
``ipp::cc``
  The ippcc library
``ipp::ch``
  The ippch library
``ipp::cv``
  The ippcv library
``ipp::dc``
  The ippdc library
``ipp::i``
  The ippi library
``ipp::s``
  The ipps library
``ipp::vm``
  The ippvm library

Result Variables
^^^^^^^^^^^^^^^^

This will define the following variables:

``IPP_DLLS``
  List of all IPP DLLs

Cache Variables
^^^^^^^^^^^^^^^

The following cache variables may also be set:

``IPP_INCLUDE_DIR``
  Include directories needed to use IPP.
``IPP_LIB_CC``
  Path to ippcc library.
``IPP_LIB_CH``
  Path to ippch library.
``IPP_LIB_CORE``
  Path to ippcore library.
``IPP_LIB_CV``
  Path to ippcv library.
``IPP_LIB_DC``
  Path to ippdc library.
``IPP_LIB_I``
  Path to ippi library.
``IPP_LIB_S``
  Path to ipps library.
``IPP_LIB_VM``
  Path to ippvm library.

#]=======================================================================]

find_path(IPP_INCLUDE_DIR ipp.h
    PATHS ${IPP_ROOT}/include)

if(UNIX)
    SET(CMAKE_FIND_LIBRARY_PREFIXES "lib")
    SET(CMAKE_FIND_LIBRARY_SUFFIXES ".so")
    set(IPP_LIB_DIR ${IPP_ROOT}/lib/intel64)
else()
    SET(CMAKE_FIND_LIBRARY_PREFIXES "")
    SET(CMAKE_FIND_LIBRARY_SUFFIXES ".lib" ".dll")
    set(IPP_LIB_DIR ${IPP_ROOT}/lib/intel64)
    set(IPP_DLL_DIR ${IPP_ROOT}/redist/intel64)
    file(GLOB IPP_DLLS "${IPP_DLL_DIR}/*.dll")
endif()

function(find_ipp_library ipp_component)
    string(TOUPPER ${ipp_component} ipp_component_uppercase)
    set(IPP_LIB_NAME IPP_LIB_${ipp_component_uppercase})
    find_library(${IPP_LIB_NAME} ipp${ipp_component} ${IPP_LIB_DIR})
    
    if(${IPP_LIB_NAME} AND NOT TARGET ipp::${ipp_component})
        add_library(ipp::${ipp_component} UNKNOWN IMPORTED)
        set_target_properties(ipp::${ipp_component} PROPERTIES
            IMPORTED_LOCATION "${${IPP_LIB_NAME}}"
            INTERFACE_INCLUDE_DIRECTORIES "${IPP_INCLUDE_DIR}"
        )
    endif()
endfunction()

foreach(component core cc ch cv dc i s vm)
    find_ipp_library(${component})
endforeach()

