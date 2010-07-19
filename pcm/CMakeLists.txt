set( LIBPCM_F90 ep_com.f90 ewaldf.f90 )
set( LIBPCM_HEADERS clj.h ewald.h lennard-jones.h functional.h )
set( LIBPCM_SOURCES clj.cc ewald.cc lennard-jones.cc erfc.cc functional.cc )

add_library (pcm SHARED ${LIBPCM_SOURCES} ${LIBPCM_HEADERS} ${LIBPCM_F90})
target_link_libraries(pcm crystal)
if(MPI_FOUND)
  include_directories(${MPI_INCLUDE_PATH})
endif(MPI_FOUND)

if(NOT noinstall_headers)
  install(FILES ${LIBPCM_HEADERS} DESTINATION include/LaDa/pcm)
endif(NOT noinstall_headers)
install(TARGETS pcm DESTINATION lib64)

if(python_bindings)
  add_subdirectory(python)
endif(python_bindings)
