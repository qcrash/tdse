#!/bin/bash

if ! id -u docker >/dev/null 2>&1; then
    useradd -m docker
fi

if [[ -n $HOST_UID ]]; then
    usermod -u $HOST_UID docker >/dev/null
fi
if [[ -n $HOST_GID ]]; then
    groupmod -g $HOST_GID docker >/dev/null
fi

source premake.sh

exec gosu docker "$@"
