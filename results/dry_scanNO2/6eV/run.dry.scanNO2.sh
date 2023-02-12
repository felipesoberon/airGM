#!/bin/bash

Te=6.0

for NO2 in 200 100 50 0
do
    
    airGM2.1 -totaltime 1E-8 -dt  2.5E-15 -Te $Te -[NO2] $NO2 -[H2O] 0
    airGM2.1 -totaltime 1E-5 -dt  2.5E-12
    airGM2.1 -totaltime 1E-4 -dt    5E-12
    airGM2.1 -totaltime 1E-3 -dt   10E-12
    airGM2.1 -totaltime 1E-2 -dt   20E-12
    airGM2.1 -totaltime 1E-1 -dt 0.25E-9
    airGM2.1 -totaltime 1    -dt  0.5E-9
    airGM2.1 -totaltime 10   -dt    1E-9
    
    gnuplot plot.gnuplot
    mv densities.pdf "densities$NO2.pdf"
    mv output.csv "output$NO2.csv"
    
done 

