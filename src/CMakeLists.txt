add_subdirectory("FD_Template")
add_subdirectory("minres-ijb")


add_library(ewald-minres STATIC fft.f90 charge_grid_module.f90 vec_utilities.f90 sort.f90 lattice.f90 ffp_module.f90 real_space_module.f90 grid_io_module.f90 ssp_module.f90 domains_module.f90 )
add_dependencies(ewald-minres minres-ijb fd)

target_include_directories(ewald-minres PUBLIC ${CMAKE_BINARY_DIR}/src/minres-ijb)
target_include_directories(ewald-minres PUBLIC ${FiniteD_BINARY_DIR}/src)

add_executable(test-ewald test.f90)
target_link_libraries(test-ewald ewald-minres minres-ijb fd fftw3)
