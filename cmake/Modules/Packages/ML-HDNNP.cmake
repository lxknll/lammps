find_package(N2P2 QUIET)
if(N2P2_FOUND)
  set(DOWNLOAD_N2P2_DEFAULT OFF)
else()
  set(DOWNLOAD_N2P2_DEFAULT ON)
endif()
option(DOWNLOAD_N2P2 "Download n2p2 library instead of using an already installed one)" ${DOWNLOAD_N2P2_DEFAULT})

if(DOWNLOAD_N2P2)
  # 1. Point the URL to the 'main' branch tarball
  set(N2P2_URL "https://github.com/CompPhysVienna/n2p2/archive/refs/heads/main.tar.gz" CACHE STRING "URL for n2p2 tarball")
  
  # 2. Comment out or remove the SHA256 checksum so CMake doesn't fail on new commits
  # set(N2P2_SHA256 "4acaa255632a7b9811d7530fd52ac7dd0bb3a8e3a3cf8512beadd29b62c1bfef" CACHE STRING "SHA256 checksum of N2P2 tarball")
  
  mark_as_advanced(N2P2_URL)
  # mark_as_advanced(N2P2_SHA256)
  
  GetFallbackURL(N2P2_URL N2P2_FALLBACK)

  find_package(Eigen3 QUIET)

  # adjust settings from detected compiler to compiler platform in n2p2 library
  if((CMAKE_CXX_COMPILER_ID STREQUAL "Clang") OR (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang"))
    set(N2P2_COMP llvm)
    set(N2P2_CXX_STD "-std=c++11")
  elseif((CMAKE_CXX_COMPILER_ID STREQUAL "Intel") OR (CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM"))
    set(N2P2_COMP intel)
    set(N2P2_CXX_STD "-std=c++11")
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(N2P2_COMP gnu)
    set(N2P2_CXX_STD "-std=gnu++11")
  elseif((CMAKE_CXX_COMPILER_ID STREQUAL "PGI") OR (CMAKE_CXX_COMPILER_ID STREQUAL "NVHPC"))
    set(N2P2_COMP gnu)
    set(N2P2_CXX_STD "--c++11")
  else() # default
    set(N2P2_COMP "")
  endif()

  # pass on archive creator command. prefer compiler specific version, if set.
  if(CMAKE_CXX_COMPILER_AR)
    set(N2P2_AR ${CMAKE_CXX_COMPILER_AR})
  else()
    set(N2P2_AR ${CMAKE_AR})
  endif()

  # adjust compilation of n2p2 library to whether MPI is requested in LAMMPS or not
  if(NOT BUILD_MPI)
    set(N2P2_PROJECT_OPTIONS "-DN2P2_NO_MPI")
  else()
    get_target_property(N2P2_MPI_INCLUDE MPI::MPI_CXX INTERFACE_INCLUDE_DIRECTORIES)
    foreach (_INCL ${N2P2_MPI_INCLUDE})
      set(N2P2_PROJECT_OPTIONS "${N2P2_PROJECT_OPTIONS} -I${_INCL}")
    endforeach()
  endif()

  # [NEW] Inject Eigen3 include directories using your specific EasyBuild environment variable
  set(N2P2_PROJECT_OPTIONS "${N2P2_PROJECT_OPTIONS} -I$ENV{EBROOTEIGEN}/include")

  # prefer GNU make, if available. N2P2 lib seems to need it.
  find_program(N2P2_MAKE NAMES gmake make)

  string(TOUPPER "${CMAKE_BUILD_TYPE}" BTYPE)
  
  # [NEW] Dynamically set OpenMP and MKL flags based on the detected compiler
  if((CMAKE_CXX_COMPILER_ID STREQUAL "Intel") OR (CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM"))
    set(N2P2_OMP_FLAG "-qopenmp")
    set(N2P2_MKL_FLAG "-qmkl")
    set(N2P2_MKL_LDFLAGS "-qmkl")
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(N2P2_OMP_FLAG "-fopenmp")
    set(N2P2_MKL_FLAG "-m64 -I$ENV{MKLROOT}/include")
    # GCC requires explicit linking against the MKL libraries
    set(N2P2_MKL_LDFLAGS "-L$ENV{MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl")
  else()
    set(N2P2_OMP_FLAG "")
    set(N2P2_MKL_FLAG "")
    set(N2P2_MKL_LDFLAGS "")
  endif()

  set(N2P2_BUILD_FLAGS "${CMAKE_SHARED_LIBRARY_CXX_FLAGS} ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${BTYPE}} ${N2P2_CXX_STD} ${N2P2_OMP_FLAG} ${N2P2_MKL_FLAG}")
  
  set(N2P2_BUILD_OPTIONS INTERFACES=LAMMPS COMP=${N2P2_COMP} "PROJECT_OPTIONS=${N2P2_PROJECT_OPTIONS}" "PROJECT_DEBUG="
    "PROJECT_CC=${CMAKE_CXX_COMPILER}" "PROJECT_MPICC=${CMAKE_CXX_COMPILER}" "PROJECT_CFLAGS=${N2P2_BUILD_FLAGS}"
    "PROJECT_LDFLAGS=${N2P2_MKL_LDFLAGS}"
    "PROJECT_AR=${N2P2_AR}" "APP_CORE=nnp-convert" "APP_TRAIN=nnp-train" "APP=nnp-convert")
  message(STATUS "N2P2 BUILD OPTIONS: ${N2P2_BUILD_OPTIONS}")

  find_program(HAVE_SED sed)
  if(NOT HAVE_SED)
    message(FATAL_ERROR "Must have 'sed' program installed to compile 'n2p2' library for ML-HDNNP package")
  endif()

  include(ExternalProject)
  ExternalProject_Add(n2p2_build
    GIT_REPOSITORY "https://github.com/CompPhysVienna/n2p2"
    GIT_TAG "master"
    GIT_SHALLOW YES
    GIT_PROGRESS YES
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    PATCH_COMMAND sed -i -e "s/\\(MPI_\\(P\\|Unp\\)ack(\\)/\\1(void *) /" src/libnnpif/LAMMPS/InterfaceLammps.cpp
    BUILD_COMMAND ${N2P2_MAKE} -C <SOURCE_DIR>/src -f makefile libnnpif ${N2P2_BUILD_OPTIONS}
    BUILD_ALWAYS YES
    INSTALL_COMMAND ""
    BUILD_IN_SOURCE 1
    LOG_BUILD ON
    SOURCE_SUBDIR src/
    BUILD_BYPRODUCTS <SOURCE_DIR>/lib/libnnp.a <SOURCE_DIR>/lib/libnnpif.a
    )

  ExternalProject_get_property(n2p2_build SOURCE_DIR)
  add_library(LAMMPS::N2P2::LIBNNP UNKNOWN IMPORTED)
  set_target_properties(LAMMPS::N2P2::LIBNNP PROPERTIES
    IMPORTED_LOCATION "${SOURCE_DIR}/lib/libnnp.a"
    INTERFACE_INCLUDE_DIRECTORIES "${SOURCE_DIR}/include")
    
  add_library(LAMMPS::N2P2::LIBNNPIF UNKNOWN IMPORTED)
  set_target_properties(LAMMPS::N2P2::LIBNNPIF PROPERTIES
    IMPORTED_LOCATION "${SOURCE_DIR}/lib/libnnpif.a"
    INTERFACE_INCLUDE_DIRECTORIES "${SOURCE_DIR}/include")
    
  if(BUILD_MPI)
    set_target_properties(LAMMPS::N2P2::LIBNNPIF PROPERTIES
      INTERFACE_LINK_LIBRARIES MPI::MPI_CXX)
    if((CMAKE_SYSTEM_NAME STREQUAL Windows) AND CMAKE_CROSSCOMPILING)
      add_dependencies(LAMMPS::N2P2::LIBNNPIF MPI::MPI_CXX)
    endif()
  endif()

  add_library(LAMMPS::N2P2 INTERFACE IMPORTED)
  set_property(TARGET LAMMPS::N2P2 PROPERTY
    INTERFACE_LINK_LIBRARIES LAMMPS::N2P2::LIBNNPIF LAMMPS::N2P2::LIBNNP)
  target_link_libraries(lammps PRIVATE LAMMPS::N2P2)

  find_package(GSL REQUIRED)
  target_link_libraries(lammps PRIVATE GSL::gsl)

  add_dependencies(LAMMPS::N2P2 n2p2_build)
  file(MAKE_DIRECTORY "${SOURCE_DIR}/include")
  file(MAKE_DIRECTORY "${SOURCE_DIR}/lib")
else()
  find_package(N2P2)
  if(NOT N2P2_FOUND)
    message(FATAL_ERROR "n2p2 not found, help CMake to find it by setting N2P2_DIR, or set DOWNLOAD_N2P2=ON to download it")
  endif()
  target_link_libraries(lammps PRIVATE N2P2::N2P2)
  include(${N2P2_CMAKE_EXTRAS})
endif()
