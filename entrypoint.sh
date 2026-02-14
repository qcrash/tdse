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

case $TARGETARCH in
    "amd64")
	source  /opt/intel/oneapi/setvars.sh > /dev/null
	;;
    "aarch64")
	echo "This is" $TARGETARCH
	;;
esac
       
exec gosu docker "$@"
