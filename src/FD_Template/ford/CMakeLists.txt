find_package(Ford)
if (Ford_FOUND)
  add_custom_target(doc
    COMMAND ${Ford_EXECUTABLE} "${FiniteD_FORD_DIR}/fd_lib.md"
    )
endif()
