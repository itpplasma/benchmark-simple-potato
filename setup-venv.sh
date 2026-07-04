#!/bin/bash
# Create and populate a benchmark-local Python virtual environment.
# Usage: ./setup-venv.sh [--recreate] [--simple-dir PATH] [--neort-dir PATH]
#
# Options:
#   --recreate         Delete an existing .venv before creating a fresh one.
#   --simple-dir PATH  Path to the SIMPLE checkout (default: ../SIMPLE).
#   --neort-dir PATH   Path to the NEO-RT checkout (default: ../NEO-RT).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
VENV_DIR="${SCRIPT_DIR}/.venv"
SIMPLE_DIR="${SCRIPT_DIR}/../SIMPLE"
NEORT_DIR="${SCRIPT_DIR}/../NEO-RT"
RECREATE=0

while [ "$#" -gt 0 ]; do
    case "$1" in
        --recreate)
            RECREATE=1
            shift
            ;;
        --simple-dir)
            if [ "$#" -lt 2 ]; then
                echo "--simple-dir requires a path argument" >&2
                exit 1
            fi
            SIMPLE_DIR="$2"
            shift 2
            ;;
        --neort-dir)
            if [ "$#" -lt 2 ]; then
                echo "--neort-dir requires a path argument" >&2
                exit 1
            fi
            NEORT_DIR="$2"
            shift 2
            ;;
        -h|--help)
            sed -n '2,/^$/s/^# //p' "$0"
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            exit 1
            ;;
    esac
done

if [ "$RECREATE" -eq 1 ] && [ -d "$VENV_DIR" ]; then
    echo "Removing existing virtual environment ${VENV_DIR} ..."
    rm -rf "$VENV_DIR"
fi

if [ -d "$VENV_DIR" ]; then
    echo "Reusing virtual environment in ${VENV_DIR} ..."
else
    echo "Creating virtual environment in ${VENV_DIR} ..."
    python3 -m venv "$VENV_DIR"
fi

# shellcheck disable=SC1091
source "${VENV_DIR}/bin/activate"

echo "Upgrading pip ..."
python -m pip install --upgrade pip

echo "Installing benchmark Python dependencies ..."
python -m pip install --prefer-binary -r "${SCRIPT_DIR}/requirements.txt"

if [ ! -d "${SIMPLE_DIR}" ]; then
    echo "SIMPLE checkout not found at ${SIMPLE_DIR}. Use --simple-dir to point to it." >&2
    exit 1
fi

if [ ! -d "${NEORT_DIR}" ]; then
    echo "NEO-RT checkout not found at ${NEORT_DIR}. Use --neort-dir to point to it." >&2
    exit 1
fi

if [ -f "${SIMPLE_DIR}/requirements.txt" ]; then
    echo "Installing SIMPLE Python dependencies from ${SIMPLE_DIR} ..."
    python -m pip install --prefer-binary -r "${SIMPLE_DIR}/requirements.txt"
fi

echo "Installing SIMPLE pysimple package from ${SIMPLE_DIR} ..."
python -m pip install --no-build-isolation -e "${SIMPLE_DIR}"

if [ -f "${NEORT_DIR}/requirements.txt" ]; then
    echo "Installing NEO-RT Python dependencies from ${NEORT_DIR} ..."
    python -m pip install --prefer-binary -r "${NEORT_DIR}/requirements.txt"
else
    echo "NEO-RT requirements.txt not found at ${NEORT_DIR}." >&2
    exit 1
fi

echo ""
echo "Done. Activate with:"
echo "  source .venv/bin/activate"
