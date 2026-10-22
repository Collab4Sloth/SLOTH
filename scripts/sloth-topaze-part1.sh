#!/bin/bash
# Part 1 (machine with Internet access): gather everything needed to build
# SLOTH on Topaze (no Internet access there) into archive.tar.gz.
export ROOT_DIR=$PWD
echo "This script uses $PWD as root directory."
mkdir -p sloth-topaze-dir && cd sloth-topaze-dir
export WORK_DIR=$ROOT_DIR/sloth-topaze-dir
export MY_LOG=pratraph
export DEST_DIR=/ccc/work/cont002/den/pratraph
echo "Create the work directory: $WORK_DIR"
echo "Your login for Topaze: $MY_LOG"
echo "Your destination directory on Topaze: $DEST_DIR"

# Same MFEM spec in part 1 (mirror) and part 2 (installation)
# Only what SLOTH uses (MUMPS: optional direct solver). For more than ~2 billion
# unknowns, add "^hypre+mixedint" (64-bit global indices) in both parts.
MFEM_SPEC="mfem@4.9.0+mpi+mumps"
SLOTH_TAG=v2.2.0

cd $WORK_DIR

echo "Getting Spack ..."
if [ ! -d "spack" ] ; then
    git clone https://github.com/spack/spack.git
fi
# Since Spack 1.0 the package recipes live in a separate repository: it must
# be in the archive, Topaze cannot download it
echo "Getting Spack packages ..."
if [ ! -d "spack-packages" ] ; then
    git clone https://github.com/spack/spack-packages.git
fi
export SPACK_ROOT=$PWD/spack
echo "Spack is cloned here: $SPACK_ROOT"
rm -rf ~/.spack
source ${SPACK_ROOT}/share/spack/setup-env.sh
spack repo set --destination $WORK_DIR/spack-packages builtin

echo "Getting PLEIADES/SLOTH ($SLOTH_TAG) ..."
if [ ! -d "sloth" ] ; then
    git clone --branch $SLOTH_TAG https://github.com/Collab4Sloth/SLOTH.git sloth
fi
echo "PLEIADES/SLOTH is cloned here: $PWD/sloth"

echo "spack bootstrap mirror --binary-packages my_bootstrap"
spack bootstrap mirror --binary-packages my_bootstrap

echo "Create Spack mirror: $PWD/mirror-mfem"
spack mirror create -d mirror-mfem -D "$MFEM_SPEC %gcc@11.2.0"

cd $ROOT_DIR
echo "Create an archive: archive.tar.gz"
tar czvf archive.tar.gz sloth-topaze-dir/
echo "Send this archive to Topaze"
scp archive.tar.gz $MY_LOG@topaze.ccc.cea.fr:$DEST_DIR/