add_subdirectory("FD_Template")

add_library(ewald-minres STATIC fft.f90 charge_grid_module.f90 vec_utilities.f90 sort.f90 lattice.f90 ffp_module.f90 real_space_module.f90 grid_io_module.f90 ssp_module.f90 domains_module.f90 quadrature_base_class.f90 quadrature_trapezium.f90 halo_base_class.f90 halo_serial.f90 halo_parallel.f90 comms_base_class.f90 comms_serial.f90 comms_parallel.f90 equation_solver_base_class.f90 equation_solver_minres.f90 equation_solver_cg.f90 Ewald_3d_module.f90 equation_solver_precon_base_class.f90 equation_solver_weighted_jacobi.f90 )
add_dependencies(ewald-minres fd)

target_include_directories(ewald-minres PUBLIC ${FiniteD_BINARY_DIR}/src)

add_executable(test-ewald test.f90)
add_executable(test-ewald_mpi test_mpi.f90)

target_link_libraries(test-ewald ewald-minres fd fftw3)
target_link_libraries(test-ewald_mpi ewald-minres fd fftw3)
