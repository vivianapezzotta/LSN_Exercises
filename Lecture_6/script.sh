#!/bin/bash

echo "Script to run the Ising simulation as a function of T"

if [ -f ene.dat ]; then
    rm ene.dat
fi

if [ -f mag.dat ]; then
    rm mag.dat
fi

if [ -f heat.dat ]; then
    rm heat.dat
fi

if [ -f chi.dat ]; then
    rm chi.dat
fi

TEMP=($(seq 0.5 0.05 2.00))

for temp in "${TEMP[@]}"
do
    echo "Simulating T=$temp"
    sed -i "1s/.*/$temp/g" input.dat
    sed -i 's/\([0-9]\),/\1./g' input.dat
    ./clean.sh
    ./Monte_Carlo_ISING_1D.exe
    tail -1 output.ene.0 | awk -v var="$temp" '{print var, $3, $4}' >> ene.dat
    sed -i 's/\([0-9]\),/\1./g' ene.dat
    tail -1 output.mag.0 | awk -v var="$temp" '{print var, $3, $4}' >> mag.dat
    sed -i 's/\([0-9]\),/\1./g' mag.dat
    tail -1 output.heat.0  | awk -v var="$temp" '{print var, $3, $4}' >> heat.dat
    sed -i 's/\([0-9]\),/\1./g' heat.dat
    tail -1 output.chi.0  | awk -v var="$temp" '{print var, $3, $4}' >> chi.dat
    sed -i 's/\([0-9]\),/\1./g' chi.dat
done
