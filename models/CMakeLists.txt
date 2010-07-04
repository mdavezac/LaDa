set( LIBPCM_SOURCES 
     ep_com.f90 ewaldf.f90
     clj.cc ewald.cc lennard-jones.cc erfc.cc functional.cc
     clj.h ewald.h lennard-jones.h functional.h )

add_library (pcm STATIC ${LIBPCM_SOURCES})
target_link_libraries(pcm crystal)

