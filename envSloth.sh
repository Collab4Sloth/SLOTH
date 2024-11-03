#=============================================
#=============================================
function Print {
    echo -e "\033[1;32m$1\033[0m"
}
#=============================================
#=============================================

if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    spack load mfem
    export HYPRE_DIR=$(spack location -i hypre)
    export MPI_DIR=$(spack location -i mpi)
    export METIS_DIR=$(spack location -i metis)

elif [[ "$OSTYPE" == "darwin"* ]]; then
    export HYPRE_DIR=$(echo $(brew --prefix hypre))
    export MPI_DIR=$(echo $(brew --prefix open-mpi))
    export METIS_DIR=$(echo $(brew --prefix metis))
    export MFEM_DIR=$(echo $(brew --prefix mfem))
else
    Print "May be an update of your OS is required. Please ensure that you are using either Darwin or Linux.".
    return
fi

#=============================================
#=============================================
DEPS_ALLOWED_OPTIONS=" --with-petsc "
BUILD_ALLOWED_OPTIONS="  --coverage, --release, --optim, --debug"
built_code="Release"
clean_build="No"
for argument; do

    arg=$(echo $argument | cut -f1 -d'=')
    value=$(echo $argument | cut -f2 -d'=')

    case "$arg" in
    --release)
        built_code="Release"
        Print "Sloth built with Release compiler options "
        ;;
    --debug)
        built_code="Debug"
        Print "Sloth built with Debug compiler options "
        ;;
    --coverage)
        built_code="Coverage"
        Print "Sloth built with Coverage compiler options "
        ;;
    --optim)
        built_code="Optim"
        Print "Sloth built with Optim compiler options "
        ;;
    --clean)
        clean_build="Yes"
        Print "Previous build will be removed if it exists "
        ;;
    *)
        Print "\nERROR with $arg in shell script options"
        return
        ;;
    esac
done

#=============================================
#=============================================
#=============================================
#=============================================
if [[ "$clean_build" == "Yes" ]]; then
    Print "Delete existing build before creating a new one..."
    rm -rf *
fi
#=============================================
Print "Create a new build..."
cmake .. ${ADDITIONAL_OPTION} -DCMAKE_BUILD_TYPE=$built_code
Print "Done!"
