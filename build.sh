module load gcc
module swap intel/19.0.5
compiler_version=`cc --version 2>/dev/null`

echo "Building for Intel"
make -f Makefile.intel clean
make -f Makefile.intel

echo "Building for GCC"
module swap PrgEnv-intel PrgEnv-gnu
make -f Makefile.gcc clean
make -f Makefile.gcc

echo "Building for Cray"
module swap PrgEnv-gnu PrgEnv-cray
make -f Makefile.cray clean
make -f Makefile.cray

