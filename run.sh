#!/bin/bash
set -e

if [ -z "$1" ]; then
        sed -e "s|output1.nc|${USER}/output1.nc|" namelist.in > namelist.tmp
else
        sed -e "s|output1.nc|${USER}/output1.nc|" "$1" > namelist.tmp
fi

mkdir -p "/tmp/${USER}"
./main.exe namelist.tmp
rm -f namelist.tmp
