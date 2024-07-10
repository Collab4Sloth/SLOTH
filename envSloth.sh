
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
BUILD_ALLOWED_OPTIONS="  --coverage, --release, --debug"
built_code="Release"
clean_build="No"
for arg in "$@"; do
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
    --clean)
        clean_build="Yes"
        Print "Previous build will be removed if it exists "
        ;;
    *)
        echo "\nERROR with $arg in shell script options"
        return
        ;;
    esac
done
#=============================================
#=============================================
if [[ "$clean_build" == "Yes" ]]; then
    echo_msg "Delete existing build before creating a new one..."
    rm -rf *
fi
#=============================================
echo_msg "Create a new build..."
cmake .. ${ADDITIONAL_OPTION} -DCMAKE_BUILD_TYPE=$built_code
echo_msg "Done!"
