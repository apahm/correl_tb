cmake_minimum_required(VERSION 3.14)

set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "Build type")

set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Build type")
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release)

# Organisation of Visual Studio solution
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set_property(GLOBAL PROPERTY AUTOGEN_TARGETS_FOLDER AutoMoc)
set_property(GLOBAL PROPERTY AUTOGEN_SOURCE_GROUP AutoMoc)
set(CMAKE_UtilityFolderName "CMakeUtilityTargets")
set_property(GLOBAL PROPERTY PREDEFINED_TARGETS_FOLDER ${CMAKE_UtilityFolderName})

set(IPP_ROOT CACHE PATH "IPP library path") 

if ( NOT IPP_ROOT )
message(FATAL_ERROR "You need to provide the correct path to all specified libraries and sources."
        "It is assumed that libraries are"
        "compartable with CMake, and `find_package(<Path>)` can correctly process configuration."     
        "`<Lib>Config.cmake` must be under that path."
        "Sources must have CMakeLists.txt"
)
endif()

project(spRLI C CXX)

add_executable(spRLI src/main.cpp)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(IPP REQUIRED)

target_compile_definitions(spRLI PRIVATE 
    _USE_MATH_DEFINES
    _UNICODE
    _NOMINMAX
    NO_OPTIX
)

target_sources(spRLI PRIVATE
    src/sampa_main.c
    src/sampa.c
    src/sampa.h
    src/spa.c
    src/spa.h
    src/bird.c
    src/bird.h
)

target_include_directories(spRLI PRIVATE  
    ${IPP_INCLUDE_DIR}
)

target_link_libraries(spRLI PRIVATE
    ipp::core
    ipp::i
    ipp::s
    ipp::vm
)

    add_custom_command(
        TARGET spRLI POST_BUILD  
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${IPP_DLLS} $<TARGET_FILE_DIR:spRLI>
    ) 