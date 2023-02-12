#!/bin/bash

for Te in 5.85 5.65
do

  airGM2.1  -totaltime   1E-8 -dt  2.5E-15  -Te $Te -[H2O] 0
  airGM2.1  -totaltime   1E-5 -dt  2.5E-12  
  airGM2.1  -totaltime   1E-4 -dt    5E-12  
  airGM2.1  -totaltime   1E-3 -dt   10E-12  
  airGM2.1  -totaltime   1E-2 -dt   20E-12  
  airGM2.1  -totaltime   1E-1 -dt 0.25E-9   
  airGM2.1  -totaltime   1    -dt  0.5E-9   
  airGM2.1  -totaltime   10   -dt    1E-9   


  gnuplot plot.gnuplot
  mv densities.pdf "densities$Te.pdf"
  mv output.csv "output$Te.csv"

done 
