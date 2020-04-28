add_library(minres-ijb SHARED halo.f90 minresModule_ijb.f90)
target_link_libraries(minres-ijb fd)
target_include_directories(minres-ijb PRIVATE ${FiniteD_BINARY_DIR}/src)
