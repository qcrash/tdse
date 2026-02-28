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
	source /usr/share/modules/init/bash
	export MODULEPATH=$MODULEPATH:/opt/arm/modulefiles
	module load arm-performance-libraries
	;;
esac
       
exec gosu docker "$@"
