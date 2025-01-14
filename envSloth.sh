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
BUILD_ALLOWED_OPTIONS="  --coverage, --release, --optim, --debug, --external"
built_code="Release"
clean_build="No"
external_values=()
use_external="OFF"

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
    --external)
        Print "Sloth built with an external package"
        
        if [[ "${$(echo $value | grep -o ',' | wc -l)}" -lt 2 ]]; then
            Print "\nError: --external must contain 3 or 4 values separated by a comma."
            Print " --external=EXT_LIBDIR,EXT_LIBNAME,EXT_SRC,EXT_TEST "
            Print " EXT_DIR : path towards the external package "
            Print " EXT_LIBNAME : path of the external dynamic library linked to SLOTH"
            Print " EXT_SRC : path of the source files of the interface between the external package and SLOTH"
            Print " EXT_TEST : path of the SLOTH tests involving the external package"
            return
        fi
        export EXT_LIBDIR=$(echo $value | cut -f1-1 -d',')
        if [[ ! -d "$EXT_LIBDIR" ]]; then
            Print "\nError:  external package does not exist. Please check the data!"
            return
        fi
        export EXT_LIBNAME=$(echo $value | cut -f2-2 -d',')

        export EXT_SRC=$(echo $value | cut -f3-3 -d',')
        if [[ ! -d "$EXT_DIR" ]]; then
            Print "\nError:  source files of the interface does not exist. Please check the data!"
            return
        fi
        EXT_TEST=""
        if [[ "${$(echo $value | grep -o ',' | wc -l)}" -eq 3 ]]; then
            export EXT_TEST=$(echo $value | cut -f4-4 -d',')

            if [[ ! -d "$EXT_TEST" ]]; then
                Print "\nError: SLOTH tests of the interface does not exist. Please check the data!"
                return
            fi
        fi

        use_external='ON'

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
cmake .. ${ADDITIONAL_OPTION} -DCMAKE_BUILD_TYPE=$built_code -DSLOTH_USE_EXTERNAL=$use_external
Print "Done!"
