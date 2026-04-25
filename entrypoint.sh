#!/bin/bash
if ! id "docker" &>/dev/null; then
    useradd docker
fi
    
if [[ -n $HOST_UID ]]; then
    usermod -u $HOST_UID docker >/dev/null
fi
if [[ -n $HOST_UID ]]; then
    groupmod -g $HOST_UID docker >/dev/null
fi

source premake.sh

exec gosu docker "$@"
