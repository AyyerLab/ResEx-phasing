file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/utils)
file(GLOB utils_src "src/*.c")
foreach(src ${utils_src})
	string(REGEX REPLACE ".*/src/(.*).c" "\\1" util ${src})
	list(APPEND c_utils ${util})
	
	add_executable(${util} ${src} $<TARGET_OBJECTS:utils.o>)
	target_link_libraries(${util} GSL::gsl GSL::gslcblas)
	target_link_libraries(${util} ${FFTWF_THREADS_LIB})
	target_link_libraries(${util} ${FFTWF_LIB})
	target_link_libraries(${util} m)
	set_target_properties(${util} PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
endforeach(src)

# Add specific object libraries for different utilities
target_sources(gen_fdens PUBLIC $<TARGET_OBJECTS:fft.o> $<TARGET_OBJECTS:volume.o>)
target_sources(band_limit PUBLIC $<TARGET_OBJECTS:fft.o>)
target_sources(boost_cont PUBLIC $<TARGET_OBJECTS:fft.o>)
target_sources(calc_fsc PUBLIC $<TARGET_OBJECTS:fft.o>)
target_sources(combine PUBLIC $<TARGET_OBJECTS:fft.o>)
target_sources(create_support PUBLIC $<TARGET_OBJECTS:fft.o>)
target_sources(gen_dens PUBLIC $<TARGET_OBJECTS:fft.o>)
target_sources(liquidize PUBLIC $<TARGET_OBJECTS:fft.o>)
target_sources(sharpen PUBLIC $<TARGET_OBJECTS:volume.o>)

# Install utilities
install(TARGETS ${c_utils} DESTINATION ResEx/utils)

install(CODE "MESSAGE(\"============================================================\")")
install(CODE "MESSAGE(\"Finished installation. Add the following to your profile file:\")")
install(CODE "MESSAGE(\"export PATH=${CMAKE_INSTALL_PREFIX}/ResEx/bin:$PATH\")")
install(CODE "MESSAGE(\"============================================================\")")
