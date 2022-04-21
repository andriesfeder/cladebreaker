#!/usr/bin/env bash

# This is a wrapper around cladebreaker for packaging the Conda recipe. It has
# been inspired by a similar wrapper in Will Rowe's DRAX pipeline
# (https://github.com/will-rowe/drax) and Robert A. Petit III's Bactopia pipeline
# (https://bactopia.github.io).

# If no user input, print usage

VERSION=0.1.0
# CONTAINER_VERSION="${VERSION%.*}.x"
# CONDA_ENV=$(which bactopia | sed 's=bin/bactopia==')
# BACTOPIA_NF="${CONDA_ENV}/share/bactopia-${CONTAINER_VERSION}"


if [[ $# == 0 ]]; then

    echo "bactopia - v${VERSION}"
fi

if [[ "$1" == "version" ]] || [[ "$1" == "--version" ]]; then
    echo "bactopia ${VERSION}"

else
    if nextflow run "${BACTOPIA_NF}/main.nf" ${OPTS} "${@:1}"; then
        # bactopia finished successfully
        if [[ "$*" == *"--cleanup_workdir"* ]] && [[ "${CAN_CLEAN_UP}" -eq 1 ]]; then
            # user asked for work dir to be cleaned up
            echo "Bactopia finished successfully! Found '--cleanup_workdir' removing '${WORK_DIR}'"
            rm -rf "${WORK_DIR}"
        fi
    else
        exit $?
    fi
fi
