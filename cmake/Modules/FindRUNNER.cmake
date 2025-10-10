### Search for a locally exposed librunner.so
include(FindPackageHandleStandardArgs)

if (DEFINED ENV{RUNNER_DIR})
  # Check if RUNNER_DIR is set manually.
  set(RUNNER_DIR "${RUNNER_DIR}")

else()
  # If not, try if directory "lib/runner/RuNNer" exists.
  get_filename_component(_fullpath "${LAMMPS_LIB_SOURCE_DIR}/runner/RuNNer" REALPATH)

  if (EXISTS ${_fullpath})
    set(RUNNER_DIR "${_fullpath}")
  endif()

endif()

# Search for the RuNNer library.
find_library(RUNNER_LIB NAMES runner HINTS "${RUNNER_DIR}/build")

find_package_handle_standard_args(RUNNER DEFAULT_MSG RUNNER_DIR RUNNER_LIB)

if(RUNNER_FOUND)
  if (NOT TARGET RUNNER::RUNNER)
    add_library(RUNNER::LIB UNKNOWN IMPORTED)
    set_target_properties(RUNNER::LIB PROPERTIES IMPORTED_LOCATION ${RUNNER_LIB})
    set(RUNNER_CMAKE_EXTRAS ${RUNNER_CMAKE_EXTRA})
  endif()
endif()

mark_as_advanced(RUNNER_DIR RUNNER_LIB)
