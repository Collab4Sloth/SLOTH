#!/usr/bin/env bash
#=============================================
#=============================================
function Print {
    echo -e "\033[1;32m$1\033[0m"
}
#=============================================
#=============================================
built_code="Release"
use_external="OFF"
local_mfem_version="No"

for argument; do

    arg=$(echo $argument | cut -f1 -d'=')
    value=$(echo $argument | cut -f2 -d'=')

    case "$arg" in
    --mfem)
        local_mfem_version="Yes"
        MFEM4SLOTH=$(echo "$value")

        # Get absolute path
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
    --external)
        Print "Sloth built with an external package"
        
        if [[ "${$(echo $value | grep -o ',' | wc -l)}" -lt 2 ]]; then
            Print "\nError: --external must contain 3 or 4 values separated by a comma."
            Print " --external=EXT_LIBDIR,EXT_LIBNAME,EXT_SRC,EXT_TEST "
            Print " EXT_LIBDIR : path towards the external package "
            Print " EXT_LIBNAME : path of the external dynamic library linked to SLOTH"
            Print " EXT_SRC : path of the source files of the interface between the external package and SLOTH"
            Print " EXT_TEST : path of the SLOTH tests involving the external package"
            
            if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
                return 1  
            else
                exit 1  
            fi
        fi
        export EXT_LIBDIR=$(echo "$value" | cut -f1 -d',')
        export EXT_LIBNAME=$(echo "$value" | cut -f2 -d',')
        export EXT_SRC=$(echo "$value" | cut -f3 -d',')
        export EXT_TEST=$(echo "$value" | cut -f4 -d',')

        # Get absolute path
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
        EXT_LIBDIR=$(to_absolute_path "$EXT_LIBDIR")
        EXT_SRC=$(to_absolute_path "$EXT_SRC")
        if [[ -n "$EXT_TEST" ]]; then
            EXT_TEST=$(to_absolute_path "$EXT_TEST")
        fi
        # Validate EXT_LIBDIR
        if [[ ! -d "$EXT_LIBDIR" ]]; then
            Print "\nError: External package directory does not exist. Please check the data!"
            
            if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
                return 1  
            else
                exit 1  
            fi
        fi

        # Validate EXT_SRC
        if [[ ! -d "$EXT_SRC" ]]; then
            Print "\nError: Source files directory does not exist. Please check the data!"
            
            if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
                return 1  
            else
                exit 1  
            fi
        fi

        # Validate EXT_TEST if present
        if [[ -n "$EXT_TEST" && ! -d "$EXT_TEST" ]]; then
            Print "\nError: Tests directory does not exist. Please check the data!"
            
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
        Print " EXT_SRC=$EXT_SRC"
        Print " EXT_TEST=${EXT_TEST:-'Not Provided'}"

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



cmake ${SCRIPT_PATH} ${ADDITIONAL_OPTION} -DCMAKE_BUILD_TYPE=$built_code -DSLOTH_USE_EXTERNAL=$use_external
Print "Done!"
