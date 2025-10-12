# Enable Fortran language support
enable_language(Fortran)

# Option to force downloading RuNNer even if found
option(DOWNLOAD_RUNNER "Force download and build of RuNNer" OFF)

# Option to build RuNNer as a shared library
option(BUILD_RUNNER_SHARED "Build RuNNer as a shared/dynamic library." OFF)

# Try to find an existing installation of RuNNer
find_package(RuNNer QUIET)

# If RuNNer is not found or a download is forced, use ExternalProject
if(NOT RuNNer_FOUND OR DOWNLOAD_RUNNER)
  message(STATUS "RuNNer not found or download forced. Building from source.")

  # Include ExternalProject module
  include(ExternalProject)

  # Add any custom CMake variables required by the RuNNer build system here.
  set(RUNNER_CMAKE_ARGS
    -DUSE_MPI=ON
    -DCMAKE_Fortran_FLAGS="-fPIC"
    -DENABLE_TESTS=OFF
  )

  # Set the library file name based on the user's choice
  if(BUILD_RUNNER_SHARED)
    set(RUNNER_LIB_NAME "libRuNNer.so")
    message(STATUS "Will build RuNNer as a shared library.")
  else()
    set(RUNNER_LIB_NAME "libRuNNer.a")
    message(STATUS "Will build RuNNer as a static library.")
  endif()

  ExternalProject_Add(runner_build
    GIT_REPOSITORY "git@gitlab.com:runner-suite/runner2.git"
    GIT_TAG "main"
    GIT_SHALLOW YES
    GIT_PROGRESS YES

    # Pass CMake arguments to RuNNer's build system
    CMAKE_ARGS
      -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
      -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
      -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
      -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
      -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/runner_install
      ${RUNNER_CMAKE_ARGS}

    # Define the build and install steps
    BUILD_COMMAND ${CMAKE_COMMAND} --build <BINARY_DIR>
    INSTALL_COMMAND ${CMAKE_COMMAND} --install <BINARY_DIR>

    # Specify the location of the built library
    BUILD_BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/runner_install/RuNNer/lib/${RUNNER_LIB_NAME}
  )

  # Get the installation directory of the external project
  ExternalProject_Get_Property(runner_build INSTALL_DIR)

  # Create an IMPORTED library target for RuNNer
  add_library(RuNNer::RuNNer UNKNOWN IMPORTED)
  set_target_properties(RuNNer::RuNNer PROPERTIES
    IMPORTED_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/runner_install/RuNNer/lib/${RUNNER_LIB_NAME}"
  )

  # Add a dependency to ensure RuNNer is built before LAMMPS
  add_dependencies(lammps runner_build)

  # Link LAMMPS to the newly built RuNNer library
  target_link_libraries(lammps PRIVATE RuNNer::RuNNer ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

else()
  # If RuNNer is found, link to the existing library
  message(STATUS "Found RuNNer: ${RuNNer_LIBRARIES}")
  target_link_libraries(lammps PRIVATE RuNNer::RuNNer ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
endif()
