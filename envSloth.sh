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
function sloth_usage() {
    cat <<'USAGE'
Usage: bash /path/to/SLOTH/envSloth.sh [options] [-- <extra CMake arguments>]
       (to be run from a build directory, not from the SLOTH root directory)

Build type (default: --release):
  --release          Release, with -march=native
  --portable         Release without -march=native (generic x86-64 target on x86)
  --debug | --coverage | --optim | --minsizerel | --relwithdebinfo

Dependencies:
  --mfem=DIR         local MFEM build (instead of spack on Linux / Homebrew on macOS)
  --exprtk=DIR       ExprTk (directory containing exprtk.hpp)
  --voro=DIR         Voro++ installation prefix
  --libtorch=DIR     LibTorch installation prefix
  --external=EXT_LIBDIR,EXT_LIBNAME
                     link Sloth with an external library

Other options:
  --shared           build Sloth as a shared library
  --install=DIR      installation directory (default: ../SlothInstall)
  --np=N             number of parallel compilation jobs (default: 4)
  -h | --help        print this help

Everything after "--" is passed unchanged to CMake, e.g.:
  bash ../envSloth.sh --exprtk=../../exprtk -- -DMATOOLS_VERBOSITY_LEVEL_1=ON
USAGE
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
ADDITIONAL_ARGS=()
np="4"

# Options of this script, then (after "--") extra arguments passed to CMake
SLOTH_ARGS=()
CMAKE_EXTRA_ARGS=()
_sloth_after_dashdash="No"
for argument in "$@"; do
    if [[ "$_sloth_after_dashdash" == "Yes" ]]; then
        CMAKE_EXTRA_ARGS+=("$argument")
    elif [[ "$argument" == "--" ]]; then
        _sloth_after_dashdash="Yes"
    else
        SLOTH_ARGS+=("$argument")
    fi
done

for argument in "${SLOTH_ARGS[@]}"; do

    # "--name=value": the value is everything after the first "="
    # (paths may contain "="); "--name" alone: empty value
    arg="${argument%%=*}"
    if [[ "$argument" == *=* ]]; then
        value="${argument#*=}"
    else
        value=""
    fi

    case "$arg" in
    -h|--help)
        sloth_usage
        return 0 2>/dev/null || exit 0  # return if sourced, exit if executed
        ;;
    --np)
        if echo "$value" | grep -qE '^[0-9]+$'; then
            np="$value"
            echo "Compilation will be done in parallel with ${np} CPU"
        else
            echo "Invalid value for --np: must be an integer"
            return 1 2>/dev/null || exit 1  # return if sourced, exit if executed
        fi
        
        ;;
    --mfem)
        local_mfem_version="Yes"
        MFEM4SLOTH=$(echo "$value")

        MFEM4SLOTH=$(to_absolute_path "$MFEM4SLOTH")
        if [[ ! -d "$MFEM4SLOTH" ]]; then
            Print "\nError: "$MFEM4SLOTH" directory does not exist. Please check the path of the local MFEM version!"
            return 1 2>/dev/null || exit 1  # return if sourced, exit if executed
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
            return 1 2>/dev/null || exit 1  # return if sourced, exit if executed
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
            return 1 2>/dev/null || exit 1  # return if sourced, exit if executed
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
            return 1 2>/dev/null || exit 1  # return if sourced, exit if executed
        else
            export ADDITIONAL_OPTION=" -DCMAKE_PREFIX_PATH=$LIBTORCH "$ADDITIONAL_OPTION
            ADDITIONAL_ARGS=("-DCMAKE_PREFIX_PATH=$LIBTORCH" "${ADDITIONAL_ARGS[@]}")
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

            return 1 2>/dev/null || exit 1  # return if sourced, exit if executed
        fi
        export EXT_LIBDIR=$(echo "$value" | cut -f1 -d',')
        export EXT_LIBNAME=$(echo "$value" | cut -f2 -d',')

        EXT_LIBDIR=$(to_absolute_path "$EXT_LIBDIR")

        # Validate EXT_LIBDIR
        if [[ ! -d "$EXT_LIBDIR" ]]; then
            Print "\nError: External package directory does not exist. Please check the data!"

            return 1 2>/dev/null || exit 1  # return if sourced, exit if executed
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
            return 1 2>/dev/null || exit 1  # return if sourced, exit if executed
        fi
        export USER_INSTALL_PATH=$(to_absolute_path "$value")
        Print "Sloth wil be installed in ${USER_INSTALL_PATH}"

        ;;
    *)
        Print "\nERROR with $arg in shell script options"

        return 1 2>/dev/null || exit 1  # return if sourced, exit if executed
        ;;
    esac
done

#=============================================
#  Environment variables
#=============================================
if [[ "$local_mfem_version" == "Yes" ]]; then
    export MFEM_DIR="$MFEM4SLOTH/mfem/INSTALLDIR/"
    export HYPRE_DIR="$MFEM4SLOTH/hypre/src/hypre/"
else
    #=============================================
    #=============================================
    #  Linux and Spack
    #=============================================
    if [[ "$OSTYPE" == "linux-gnu"* ]]; then
        # Not finding spack or a package is not an error by itself: MFEM may
        # also be found by CMake through CMAKE_PREFIX_PATH (spack environment,
        # "module load", ...). Clear warnings are printed instead.
        if command -v spack &>/dev/null; then
            if ! spack load mfem 2>/dev/null; then
                Print "Note: 'spack load mfem' not possible here (no spack shell support in this shell, or mfem not found); paths from 'spack location -i' are used."
            fi
            for _sloth_pkg in MFEM_DIR:mfem HYPRE_DIR:hypre MPI_DIR:mpi; do    
                _sloth_var="${_sloth_pkg%%:*}"
                _sloth_name="${_sloth_pkg#*:}"
                if _sloth_loc=$(spack location -i "$_sloth_name" 2>/dev/null) && [[ -n "$_sloth_loc" ]]; then
                    export "$_sloth_var=$_sloth_loc"
                else
                    export "$_sloth_var="
                    Print "Warning: 'spack location -i $_sloth_name' failed (not installed, or several installations: activate a spack environment or use 'spack load <hash>'). $_sloth_var is left empty."
                fi
            done
        else
            export MPI_DIR= MFEM_DIR= HYPRE_DIR=
            Print "Warning: spack command not found. MFEM must then be found through CMAKE_PREFIX_PATH (spack environment, module load, ...) or given with --mfem=DIR."
        fi
        #=============================================
        #  Mac and Homebrew
        #=============================================
    elif [[ "$OSTYPE" == "darwin"* ]]; then
        if command -v brew &>/dev/null; then
            export HYPRE_DIR=$(brew --prefix hypre)
            export MFEM_DIR=$(brew --prefix mfem)
        else
            export HYPRE_DIR= MFEM_DIR=
            Print "Warning: brew command not found. MFEM must then be found through CMAKE_PREFIX_PATH or given with --mfem=DIR."
        fi
    else
        Print "May be an update of your OS is required. Please ensure that you are using either Darwin or Linux.".
        return 1 2>/dev/null || exit 1  # return if sourced, exit if executed
    fi
#=============================================
#=============================================
fi
#=============================================
Print "Create a new build..."

# Directory of this script (BASH_SOURCE: also correct when the script is
# sourced from bash; zsh sets $0 to the script path when sourcing)
SCRIPT_PATH=$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)
CURRENT_PATH=$(pwd)

if [ "$SCRIPT_PATH" = "$CURRENT_PATH" ] || [ "$CURRENT_PATH" = "/" ]; then
    echo "Error: Do not run the build from the root directory."
    echo "Please create a new directory to build SLOTH."
    return 1 2>/dev/null || exit 1  # return if sourced, exit if executed
else
    echo "SLOTH root path: $SCRIPT_PATH"
    echo "Current build path: $CURRENT_PATH"
fi

Print "Run CMake configuration..."
if ! cmake "${SCRIPT_PATH}" "${ADDITIONAL_ARGS[@]}" -DCMAKE_BUILD_TYPE=$built_code -DSLOTH_USE_SHARED=$use_shared -DSLOTH_USE_EXTERNAL=$use_external -DSLOTH_USE_LIBTORCH=$use_libtorch -DSLOTH_USE_EXPRTK=$use_exprtk -DSLOTH_USE_VORO=$use_voro "${CMAKE_EXTRA_ARGS[@]}"; then
    Print "\nError: CMake configuration failed (see the messages above). SLOTH is not compiled."
    return 1 2>/dev/null || exit 1  # return if sourced, exit if executed
fi

Print "Compile SLOTH..."

if ! make -j ${np}; then
    Print "\nError: compilation failed (see the messages above)."
    return 1 2>/dev/null || exit 1  # return if sourced, exit if executed
fi

Print "Done!"
