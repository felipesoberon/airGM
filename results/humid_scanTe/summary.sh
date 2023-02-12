#!/bin/bash

filename=summary.csv

rm $filename *~

head -n 1 output3.0.csv | tr -d '\n' >> $filename
echo ",Te" >> $filename


for Te in 3.0 4.0 5.0 5.25 5.50 5.75 6.0
do
    tail -n 1 output$Te.csv | tr -d '\n' >> $filename
    echo ","$Te >> $filename
done
