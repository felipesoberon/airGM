set terminal pdf
set output "summary.pdf"

set logscale y

set xlabel "Electron Temperature (eV)"
set ylabel "Density (m-3)"

set yrange [1E0:1E25]
set grid
show grid

set datafile separator ","


set title "Charged Species (e and ions)"
set key left top
plot 'summary.csv' u 56:($18+$19+$20+$21+$22+$23+$24+$25+$26+$27) \
                     title 'Ions-' w lp, \
     'summary.csv' u 56:($1+$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16) \
                     title 'Ions+' w lp, \
     'summary.csv' u 56:17 title 'e' w lp, 


#POSITIVE IONS +++++++++++++++++++++++++++++++++++

set title "Positive Nitrogen Ion Species"
set key left top
plot 'summary.csv' u 56:1 title 'N+' w lp, \
     'summary.csv' u 56:2 title 'N_2+' w lp, \
     'summary.csv' u 56:3 title 'N_3+' w lp, \
     'summary.csv' u 56:4 title 'N_4+' w lp 

set title "Positive Oxygen Ion Species"
set key left top
plot 'summary.csv' u 56:5 title 'O+' w lp, \
     'summary.csv' u 56:6 title 'O_2+' w lp, \
     'summary.csv' u 56:7 title 'O_4+' w lp 

set title "Positive Nitrogen-Oxide Ion Species"
set key left top
plot 'summary.csv' u 56:8  title 'NO+' w lp, \
     'summary.csv' u 56:9  title 'N_2O+' w lp, \
     'summary.csv' u 56:10 title 'NO_2+' w lp

set title "Positive Hydrogen Ion Species"
set key left top
plot 'summary.csv' u 56:11 title 'H+' w lp, \
     'summary.csv' u 56:12 title 'H_2+' w lp, \
     'summary.csv' u 56:13 title 'H_3+' w lp, \
     'summary.csv' u 56:14 title 'OH+' w lp, \
     'summary.csv' u 56:15 title 'H_2O+' w lp, \
     'summary.csv' u 56:16 title 'H_3O+' w lp


#NEGATIVE IONS ------------------------------------

set title "Negative Oxygen Ion Species"
set key left top
plot 'summary.csv' u 56:18 title 'O-' w lp, \
     'summary.csv' u 56:19 title 'O_2-' w lp, \
     'summary.csv' u 56:20 title 'O_3-' w lp, \
     'summary.csv' u 56:21 title 'O_4-' w lp 


set title "Negative Nitrogen-Oxide Ion Species"
set key left top
plot 'summary.csv' u 56:22 title 'NO-' w lp, \
     'summary.csv' u 56:23 title 'N_2O-' w lp, \
     'summary.csv' u 56:24 title 'NO_2-' w lp, \
     'summary.csv' u 56:25 title 'NO_3-' w lp 


set title "Negative Hydrogen Ion Species"
set key left top
plot 'summary.csv' u 56:26 title 'H-' w lp, \
     'summary.csv' u 56:27 title 'OH-' w lp, \



#NEUTRAL SPECIES ===================================

set title "Nitrogen Excited Species"
set key left top
plot 'summary.csv' u 56:28 title 'N(2D)' w lp, \
     'summary.csv' u 56:29 title 'N_2(A3Sigma)' w lp, \
     'summary.csv' u 56:30 title 'N_2(B3Pi)' w lp


set title "Atomic and Excited/Reactive Oxygen Species"
set key left top
plot 'summary.csv' u 56:31 title 'O(1D)' w lp, \
     'summary.csv' u 56:32 title 'H' w lp, \
     'summary.csv' u 56:33 title 'N' w lp, \
     'summary.csv' u 56:34 title 'O' w lp, \
     'summary.csv' u 56:35 title 'O_2(a1Delta)' w lp, \
     'summary.csv' u 56:36 title 'O_3' w lp



set title "Nitrogen-Oxide Species"
set key left top
plot 'summary.csv' u 56:37 title 'NO' w lp, \
     'summary.csv' u 56:38 title 'N_2O' w lp, \
     'summary.csv' u 56:39 title 'NO_2' w lp, \
     'summary.csv' u 56:40 title 'NO_3' w lp, \
     'summary.csv' u 56:41 title 'N_2O_3' w lp, \
     'summary.csv' u 56:42 title 'N_2O_4' w lp, \
     'summary.csv' u 56:43 title 'N_2O_5' w lp


set title "Hydrogen containing neutral Species"
set key left top
plot 'summary.csv' u 56:44 title 'H_2' w lp, \
     'summary.csv' u 56:45 title 'OH' w lp, \
     'summary.csv' u 56:46 title 'HO_2' w lp, \
     'summary.csv' u 56:47 title 'H_2O_2' w lp, \
     'summary.csv' u 56:48 title 'HNO' w lp, \
     'summary.csv' u 56:49 title 'HNO_2' w lp, \
     'summary.csv' u 56:50 title 'HNO_3' w lp


set yrange [1E19:1E26]
set title "N_2, O_2 and H_2O"
plot 'summary.csv' u 56:51 title 'N_2' w lp, \
     'summary.csv' u 56:52 title 'O_2' w lp, \
     'summary.csv' u 56:53 title 'H_2O' w lp


#PLOTS IN PPB ========================

Mdensity = 2.4E25
m3PPB = 1E9/Mdensity

set ylabel "Density (PPB)"
set yrange [1E-12:1E8]
set format y "%.1E"

set title "Reactive Species"
set key left top
plot 'summary.csv' u 56:($34*m3PPB) title 'O' w lp, \
     'summary.csv' u 56:($36*m3PPB) title 'O_3' w lp, \
     'summary.csv' u 56:($39*m3PPB) title 'NO_2' w lp, \
     'summary.csv' u 56:($40*m3PPB) title 'NO_3' w lp
     
