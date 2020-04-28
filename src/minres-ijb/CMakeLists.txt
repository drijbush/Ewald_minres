add_library(minres-ijb STATIC halo.f90 minresModule_ijb.f90)
add_dependencies(minres-ijb fd)
target_include_directories(minres-ijb PUBLIC ${FiniteD_BINARY_DIR}/src)
