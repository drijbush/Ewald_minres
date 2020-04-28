# Set up library
add_library(fd constants.f90 fd_template.f90 grid_vectors.f90 fd_laplacian_3d.f90 )

add_executable(test_fd test.f90)
target_link_libraries(test_fd fd m)
