find_program(Ford_EXECUTABLE ford QUIET DOC "Automatic documentation generator for modern Fortran programs" )

if( Ford_EXECUTABLE )

  execute_process(
    COMMAND ${Ford_EXECUTABLE} --version
    OUTPUT_VARIABLE Ford_VERSION
    ERROR_VARIABLE Ford_VERSION_BAK
    )

  if (Ford_VERSION)
    string(REGEX MATCH "[1-9]+(\.[0-9]+)+" Ford_VERSION ${Ford_VERSION})
  elseif(Ford_VERSION_BAK)
    set(Ford_VERSION ${Ford_VERSION_BAK})
    string(REGEX MATCH "[1-9]+(\.[0-9]+)+" Ford_VERSION ${Ford_VERSION})
  endif()

  find_package_handle_standard_args(
    Ford
    REQUIRED_VARS Ford_EXECUTABLE
    VERSION_VAR Ford_VERSION
    )

  mark_as_advanced(Ford_EXECUTABLE)
  mark_as_advanced(Ford_VERSION)
else()
  set( Ford_FOUND FALSE )
endif()
