#!/bin/sh

exec 3>&1 4>&2 >install.txt 2>&1

if ! command -v sage &> /dev/null
then
    if ! command -v micromamba &> /dev/null
    then
        echo "installing micromamba"
        curl micro.mamba.pm/install.sh | bash
    else
        echo "micromamba installed!"
    fi
    eval "$(micromamba shell hook --shell=bash)"
    micromamba activate sage_env
    if [ $? -ne 0 ]; then
        micromamba create -y -n sage_env sage python=3.8
        micromamba activate sage_env
    fi
fi

exec 1>&3 2>&4