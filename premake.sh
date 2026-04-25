#!/bin/bash

case $TARGETARCH in
    "amd64")
	source  /opt/intel/oneapi/setvars.sh > /dev/null
	;;
    "arm64")
	source /usr/share/modules/init/bash
	export MODULEPATH=$MODULEPATH:/opt/arm/modulefiles
	module load arm-performance-libraries
	;;
    *)
	echo "TARGETARCH not recognized"
	exit 1
	;;
esac

