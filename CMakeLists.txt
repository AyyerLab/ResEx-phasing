cmake_minimum_required(VERSION 3.6)
# Version 3.6 required for list(FILTER)

set(CMAKE_INSTALL_MESSAGE LAZY)

project(ResEx-phasing LANGUAGES C)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
	set(CMAKE_INSTALL_PREFIX $ENV{HOME}/.local CACHE PATH "Default installation path" FORCE)
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -Wall")
if(CMAKE_COMPILER_IS_GNUCC)
	if(CMAKE_C_COMPILER_VERSION VERSION_GREATER 4.4)
		#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unused-result")
		set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unused-result -Wno-format-overflow")
	endif()
endif(CMAKE_COMPILER_IS_GNUCC)

find_package(FFTW REQUIRED)
include_directories(${FFTW_INCLUDES})
find_package(OpenMP REQUIRED)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})

# Compile object files
file(GLOB OBJSOURCE "src/*.c")
list(FILTER OBJSOURCE EXCLUDE REGEX "recon.c")
foreach(src ${OBJSOURCE})
	string(REGEX REPLACE ".*/src/(.*).c" "\\1.o" obj ${src})
	add_library(${obj} OBJECT ${src})
	list(APPEND OBJECTS $<TARGET_OBJECTS:${obj}>)
endforeach(src)

# Compile recon
add_executable(recon "src/recon.c" "src/algorithm.h" ${OBJECTS})
target_link_libraries(recon GSL::gsl GSL::gslcblas)
target_link_libraries(recon fftw3f_threads)
target_link_libraries(recon fftw3f)
target_link_libraries(recon m)
set_target_properties(recon PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)

# Install recon and GUI
install(TARGETS recon DESTINATION ResEx/bin)
install(PROGRAMS ResEx DESTINATION ResEx/bin)
install(DIRECTORY gui DESTINATION ResEx)

# Compile and install utilities
add_subdirectory(utils)

