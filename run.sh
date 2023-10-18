#!/bin/bash

sed -e "s|output1.nc|${USER}/output1.nc|" namelist.in > namelist.tmp

mkdir /tmp/${USER}
./main.exe namelist.tmp
rm namelist.tmp 



