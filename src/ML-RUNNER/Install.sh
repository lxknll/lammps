# Install/unInstall package files in LAMMPS
# mode = 0/1/2 for uninstall/install/update

mode=$1

if (test $mode = 1 || test $mode = 2 || test $mode = 3) then
  echo "The ML-RUNNER package does not support the legacy build system. Please build LAMMPS with CMake instead."
  exit 1
