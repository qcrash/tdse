#!/bin/bash

if [[ -n $HOST_UID ]]; then
    usermod -u $HOST_UID docker >/dev/null
fi
if [[ -n $HOST_UID ]]; then
    groupmod -g $HOST_UID docker >/dev/null
fi

source /opt/intel/oneapi/setvars.sh > /dev/null
exec gosu docker "$@"
