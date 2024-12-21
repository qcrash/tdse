#!/bin/bash

if [[ -n $HOST_UID ]]; then
    usermod -u $HOST_UID docker >/dev/null
fi
if [[ -n $HOST_UID ]]; then
    groupmod -g $HOST_UID docker >/dev/null
fi

exec gosu docker /bin/bash --rcfile /user/docker/.bashrc -c "env"
