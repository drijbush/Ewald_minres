add_library(minres-ijb STATIC halo_base_class.f90 halo_serial.f90 minresModule_ijb.f90)
add_dependencies(minres-ijb fd)
target_include_directories(minres-ijb PUBLIC ${FiniteD_BINARY_DIR}/src)
