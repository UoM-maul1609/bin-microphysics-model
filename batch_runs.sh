#! /bin/bash
ARRAY1=(0.354875716148452   1.027155925684173   1.699436135219894   2.371716344755614 3.043996554291334   3.716276763827054) # cloud base
ELEMENTS1=${#ARRAY1[@]} # elements in first array

for (( i=0;i<$ELEMENTS1;i++)); do
    # Runs with the hm process switched on:
    sed -e "s/winit=.6,/winit=${ARRAY1[${i}]},/" namelist.in >/tmp/namelist.tmp
    cp /tmp/namelist.tmp ./namelist.run

    echo ${ARRAY1[${i}]} ' no entrainment'
    ./main.exe namelist.run > /tmp/std.out
    mv /tmp/output1.nc /tmp/output_${i}_1_bmm.nc

    sed -e "s/adiabatic_prof=.true.,/adiabatic_prof=.false.,/" /tmp/namelist.tmp > namelist.run

    echo ${ARRAY1[${i}]} ' with entrainment'
    ./main.exe namelist.run > /tmp/std.out
    mv /tmp/output1.nc /tmp/output_${i}_2_bmm.nc


done

