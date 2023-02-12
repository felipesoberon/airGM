#!/bin/bash

Te=6.0
tp=1E-9

#          time    Te(eV)  dt(s)    [H2O](m-3) tplasma
airGM2     1E-8    $Te     2.5E-15  0          $tp
airGM2     1E-5    $Te     2.5E-12  0          $tp
airGM2     1E-4    $Te     5.0E-12  0          $tp
airGM2     1E-3    $Te    10.0E-12  0          $tp
airGM2     1E-2    $Te    20.0E-12  0          $tp
airGM2     1E-1    $Te    50.0E-12  0          $tp
airGM2     1       $Te    50.0E-12  0          $tp
airGM2     10      $Te       1E-9   0          $tp
airGM2     100     $Te       1E-9   0          $tp
