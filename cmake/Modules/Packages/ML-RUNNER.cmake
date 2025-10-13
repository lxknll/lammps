# Enable Fortran language support
enable_language(Fortran)

# Build arguments.
option(DOWNLOAD_RUNNER "Force download and build of RuNNer. If this is OFF\
  the arguments RUNNER_LIB_DIR, RUNNER_LIB_NAME and RUNNER_SHARED_LIB must\
  be set to meaningful values." OFF
)
option(RUNNER_SHARED_LIB "Use pre-compiled shared/dynamic RuNNer library. Only \
  considered if DOWNLOAD_RUNNER is OFF." OFF
)
set(RUNNER_LIB_DIR "$ENV{HOME}/.local/lib" CACHE STRING "Root directory of a \
  pre-compiled RuNNer installation. Only considered if DOWNLOAD_RUNNER is OFF."
)
set(RUNNER_LIB_NAME "libRuNNer" CACHE STRING "Name of the RuNNer library \
  (excluding file extension, this is controlled by RUNNER_SHARED_LIB). Only \
  considered if DOWNLOAD_RUNNER is OFF."
)

if(DOWNLOAD_RUNNER)
  message(STATUS "DOWNLOAD_RUNNER is ON. Building RuNNer from source as a static library.")

  # Include ExternalProject module
  include(ExternalProject)

  # Add any custom CMake variables required by the RuNNer build system here.
  set(RUNNER_CMAKE_ARGS
    -DUSE_MPI=ON
    -DCMAKE_Fortran_FLAGS="-fPIC"
    -DENABLE_TESTS=no
  )

  # Force using the static library. RuNNer's cmake build always produces libRuNNer.a.
  set(RUNNER_LIB_NAME "libRuNNer.a")

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

  # Create an IMPORTED library target for RuNNer
  add_library(RuNNer::RuNNer STATIC IMPORTED)
  set_target_properties(RuNNer::RuNNer PROPERTIES
    IMPORTED_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/runner_install/RuNNer/lib/${RUNNER_LIB_NAME}"
  )

  # Add a dependency to ensure RuNNer is built before the main target
  add_dependencies(lammps runner_build)

else()
  message(STATUS "DOWNLOAD_RUNNER is OFF. Looking for a pre-compiled RuNNer.")

  # Check if the directory specified by RUNNER_LIB_DIR exists
  if(NOT IS_DIRECTORY ${RUNNER_LIB_DIR})
    message(FATAL_ERROR "The directory specified by RUNNER_LIB_DIR does not exist: ${RUNNER_LIB_DIR}")
  endif()

  # Determine the correct file extension based on whether the library is shared or not.
  if(RUNNER_SHARED_LIB)
    set(RUNNER_LIB_EXT ".so")
    set(RUNNER_LIB_TYPE SHARED)
  else()
    set(RUNNER_LIB_EXT ".a")
    set(RUNNER_LIB_TYPE STATIC)
  endif()

  set(RUNNER_LIB_FULL_NAME ${RUNNER_LIB_NAME}${RUNNER_LIB_EXT})
  message(STATUS "Looking for the ${RUNNER_LIB_TYPE} library: ${RUNNER_LIB_FULL_NAME}")

  # Find the RuNNer library in the specified path
  find_library(RuNNer_LIBRARY
    NAMES ${RUNNER_LIB_FULL_NAME}
    HINTS "${RUNNER_LIB_DIR}"
    NO_DEFAULT_PATH
  )

  if(RuNNer_LIBRARY)
    message(STATUS "Found RuNNer library: ${RuNNer_LIBRARY}")

    # Create an IMPORTED library target for the found RuNNer library
    add_library(RuNNer::RuNNer ${RUNNER_LIB_TYPE} IMPORTED)
    set_target_properties(RuNNer::RuNNer PROPERTIES
      IMPORTED_LOCATION "${RuNNer_LIBRARY}"
    )
  else()
    message(FATAL_ERROR "Could not find the RuNNer library in the specified \
      RUNNER_LIB_DIR. Searched for ${RUNNER_LIB_FULL_NAME} at ${RUNNER_LIB_DIR}"
    )
  endif()

endif()

# Link lammps to the found RuNNer library
target_link_libraries(lammps PRIVATE RuNNer::RuNNer ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
