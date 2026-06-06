#!/usr/bin/env bash
#=============================================
#=============================================
function Print {
    echo -e "\033[1;32m$1\033[0m"
}

function to_absolute_path() {
    local path=$1
    if [[ ! "$path" =~ ^/ ]]; then
        if command -v realpath &>/dev/null; then
            echo "$(realpath "$path")"
        else
            echo "$(cd "$path" && pwd)"
        fi
    else
        echo "$path"
    fi
}
#=============================================
#=============================================
built_code="Release"
use_external="OFF"
use_libtorch="OFF"
use_shared="OFF"
use_exprtk="OFF"
use_voro="OFF"
local_mfem_version="No"
export USER_INSTALL_PATH=""
ADDITIONAL_OPTION=""
np="4"

for argument; do

    arg=$(echo $argument | cut -f1 -d'=')
    value=$(echo $argument | cut -f2 -d'=')

    case "$arg" in
    --np)
        if echo "$value" | grep -qE '^[0-9]+$'; then
            np="$value"
            echo "Compilation will be done in parallel with ${np} CPU"
        else
            echo "Invalid value for --np: must be an integer"
            exit 1
        fi
        
        ;;
    --mfem)
        local_mfem_version="Yes"
        MFEM4SLOTH=$(echo "$value")

        MFEM4SLOTH=$(to_absolute_path "$MFEM4SLOTH")
        if [[ ! -d "$MFEM4SLOTH" ]]; then
            Print "\nError: "$MFEM4SLOTH" directory does not exist. Please check the path of the local MFEM version!"
            if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
                return 1
            else
                exit 1
            fi
        else
            Print "Sloth built with a local MFEM version: $MFEM4SLOTH"
            export MFEM4SLOTH=$MFEM4SLOTH
        fi
        ;;
    --shared)
        use_shared="ON"
        Print "Shared library of Sloth will be created"
        ;;
    --release)
        built_code="Release"
        Print "Sloth built with Release compiler options [with -march=native]"
        ;;
    --portable)
        built_code="PortableRelease"
        Print "Sloth built with PortableRelease compiler options [Release with -march=x86_64]"
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
    --minsizerel)
        built_code="MinSizeRel"
        Print "Sloth built with MinSizeRel compiler options "
        ;;
    --relwithdebinfo)
        built_code="RelWithDebInfo"
        Print "Sloth built with RelWithDebInfo compiler options "
        ;;
    --voro)
        export VORO=$(to_absolute_path "$value")
        if [[ ! -d "$VORO" ]]; then
            Print "\nError: "$VORO" directory does not exist. Please check the path of the local voro++ version!"
            if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
                return 1
            else
                exit 1
            fi
        else
            Print "Sloth built with voro++ "
            # Set voro++ flag
            use_voro='ON'
        fi
        ;;
    --exprtk)
        # export ExprTk_INCLUDE_DIR=$(spack location -i exprtk)/include/exprtk

        export ExprTk_INCLUDE_DIR=$(to_absolute_path "$value")
        if [[ ! -d "$ExprTk_INCLUDE_DIR" ]]; then
            Print "\nError: "$ExprTk_INCLUDE_DIR" directory does not exist. Please check the path of the local ExprTk version!"
            if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
                return 1
            else
                exit 1
            fi
        else
            Print "Sloth built with ExprTk "
            # Set libtorch flag
            use_exprtk='ON'
        fi

        ;;
    --libtorch)
        export LIBTORCH=$(to_absolute_path "$value")
        if [[ ! -d "$LIBTORCH" ]]; then
            Print "\nError: "$LIBTORCH" directory does not exist. Please check the path of the local LIBTORCH version!"
            if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
                return 1
            else
                exit 1
            fi
        else
            export ADDITIONAL_OPTION=" -DCMAKE_PREFIX_PATH=$LIBTORCH "$ADDITIONAL_OPTION
            Print "Sloth built with libTorch "
            # Set libtorch flag
            use_libtorch='ON'
        fi

        ;;
    --external)
        Print "Sloth built with an external package that requires linking to an external library."
        count=$(echo "$value" | grep -o ',' | wc -l)

        if [[ "${count}" -lt 1 ]]; then
            Print "\nError: --external must contain 2 values separated by a comma."
            Print " --external=EXT_LIBDIR,EXT_LIBNAME,EXT_SRC,EXT_TEST "
            Print " EXT_LIBDIR : path towards the external package "
            Print " EXT_LIBNAME : path of the external dynamic library linked to SLOTH"

            if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
                return 1
            else
                exit 1
            fi
        fi
        export EXT_LIBDIR=$(echo "$value" | cut -f1 -d',')
        export EXT_LIBNAME=$(echo "$value" | cut -f2 -d',')

        EXT_LIBDIR=$(to_absolute_path "$EXT_LIBDIR")

        # Validate EXT_LIBDIR
        if [[ ! -d "$EXT_LIBDIR" ]]; then
            Print "\nError: External package directory does not exist. Please check the data!"

            if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
                return 1
            else
                exit 1
            fi
        fi

        # Set external flag
        use_external='ON'

        Print "\nConfiguration successful:"
        Print " EXT_LIBDIR=$EXT_LIBDIR"
        Print " EXT_LIBNAME=$EXT_LIBNAME"

        ;;
    --install)
        if [ -z "$value" ] || [ "${value#--}" != "$value" ]; then
            echo "Error: --install requires a valid path"
            exit 1
        fi
        export USER_INSTALL_PATH=$(to_absolute_path "$value")
        Print "Sloth wil be installed in ${USER_INSTALL_PATH}"

        ;;
    *)
        Print "\nERROR with $arg in shell script options"

        if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
            return 1
        else
            exit 1
        fi
        ;;
    esac
done

#=============================================
#  Environment variables
#=============================================
if [[ "$local_mfem_version" == "Yes" ]]; then
    export MFEM_DIR="$MFEM4SLOTH/mfem/INSTALLDIR/"
    export HYPRE_DIR="$MFEM4SLOTH/hypre/src/hypre/"
    export METIS_DIR="$MFEM4SLOTH/metis-4.0/"
    export SuiteSparse_DIR="$MFEM4SLOTH/SuiteSparse/"
else
    #=============================================
    #=============================================
    #  Linux and Spack
    #=============================================
    if [[ "$OSTYPE" == "linux-gnu"* ]]; then
        spack load mfem
        export MFEM_DIR=$(spack location -i mfem)
        export HYPRE_DIR=$(spack location -i hypre)
        export MPI_DIR=$(spack location -i mpi)
        export METIS_DIR=$(spack location -i metis)
        #=============================================
        #  Mac and Homebrew
        #=============================================
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
fi
#=============================================
Print "Create a new build..."

SCRIPT_PATH=$(cd "$(dirname "$0")" && pwd)
CURRENT_PATH=$(pwd)

if [ "$SCRIPT_PATH" = "$CURRENT_PATH" ] || [ "$CURRENT_PATH" = "/" ]; then
    echo "Error: Do not run the build from the root directory."
    echo "Please create a new directory to build SLOTH."
    exit 1
else
    echo "SLOTH root path: $SCRIPT_PATH"
    echo "Current build path: $CURRENT_PATH"
fi

Print "Run CMake configuration..."
cmake ${SCRIPT_PATH} ${ADDITIONAL_OPTION}  -DCMAKE_BUILD_TYPE=$built_code -DSLOTH_USE_SHARED=$use_shared -DSLOTH_USE_EXTERNAL=$use_external  -DSLOTH_USE_LIBTORCH=$use_libtorch  -DSLOTH_USE_EXPRTK=$use_exprtk  -DSLOTH_USE_VORO=$use_voro

Print "Compile SLOTH..."

make -j ${np}

Print "Done!"
