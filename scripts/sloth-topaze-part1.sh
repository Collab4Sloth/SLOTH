#!/bin/sh
#PS1='\[\033[1;36m\]\u\[\033[1;31m\]@\[\033[1;32m\]\h:\[\033[1;35m\]\w\[\033[1;31m\]\$\[\033[0m\] '
export ROOT_DIR=$PWD
echo "This script uses $PWD as root directory."   
mkdir -p sloth-topaze-dir && cd sloth-topaze-dir
export WORK_DIR=$ROOT_DIR/sloth-topaze-dir
export MY_LOG=pratraph
export DEST_DIR=/ccc/work/cont002/den/pratraph
echo "Create the work directory: $WORK_DIR"
echo "Your login for Topaze: $MY_LOG"
echo "Your destination directory on Topaze: $DEST_DIR"
   

cd $WORK_DIR

echo "Getting Spack ..."
if [ ! -d "spack" ] ; then
 	git clone https://github.com/spack/spack.git
fi
export SPACK_ROOT=$PWD/spack
echo "Spack is clone here: $SPACK_ROOT"
rm -r ~/.spack
source ${SPACK_ROOT}/share/spack/setup-env.sh
echo "Getting PLEIADES/SLOTH ..."
if [ ! -d "sloth" ] ; then
	git clone https://www-git-cad.intra.cea.fr/DEC/collaboratif/ci230846/DEV_PROJECT/sloth.git
fi
echo "PLEIADES/SLOTH is clone here: $PWD/sloth"
echo "spack bootstrap mirror --binary-packages my_bootstrap"
spack bootstrap mirror --binary-packages my_bootstrap
echo "Create Spack mirror: $PWD/mirror-mfem"
spack mirror create -d mirror-mfem -D gcc@11.2.0 mfem+mpi+debug+openmp+petsc+strumpack+suite-sparse+sundials+superlu-dist%gcc@11.2.0
cd $ROOT_DIR
echo "Create an archive: archive.tar.gz"
tar cvf archive.tar.gz sloth-topaze-dir/
echo "Send this archive to Topaze"
scp archive.tar.gz $MY_LOG@topaze.ccc.cea.fr:$DEST_DIR/
