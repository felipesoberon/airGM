set terminal pdf
set output "densities.pdf"

set logscale y
set logscale x

set xlabel "Time (s)"
set ylabel "Density (m-3)"

set yrange [1E0:1E21]
set grid
show grid

set datafile separator ","


set title "Charged Species (e and ions)"
set key left top
plot 'output.csv' u 54:($18+$19+$20+$21+$22+$23+$24+$25+$26+$27) \
                     title 'Ions-' w l, \
     'output.csv' u 54:($1+$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16) \
                     title 'Ions+' w l, \
     'output.csv' u 54:17 title 'e' w l, 


#POSITIVE IONS +++++++++++++++++++++++++++++++++++

set title "Positive Nitrogen Ion Species"
set key left top
plot 'output.csv' u 54:1 title 'N+' w l, \
     'output.csv' u 54:2 title 'N_2+' w l, \
     'output.csv' u 54:3 title 'N_3+' w l, \
     'output.csv' u 54:4 title 'N_4+' w l 

set title "Positive Oxygen Ion Species"
set key left top
plot 'output.csv' u 54:5 title 'O+' w l, \
     'output.csv' u 54:6 title 'O_2+' w l, \
     'output.csv' u 54:7 title 'O_4+' w l 

set title "Positive Nitrogen-Oxide Ion Species"
set key left top
plot 'output.csv' u 54:8  title 'NO+' w l, \
     'output.csv' u 54:9  title 'N_2O+' w l, \
     'output.csv' u 54:10 title 'NO_2+' w l

set title "Positive Hydrogen Ion Species"
set key left top
plot 'output.csv' u 54:11 title 'H+' w l, \
     'output.csv' u 54:12 title 'H_2+' w l, \
     'output.csv' u 54:13 title 'H_3+' w l, \
     'output.csv' u 54:14 title 'OH+' w l, \
     'output.csv' u 54:15 title 'H_2O+' w l, \
     'output.csv' u 54:16 title 'H_3O+' w l


#NEGATIVE IONS ------------------------------------

set title "Negative Oxygen Ion Species"
set key left top
plot 'output.csv' u 54:18 title 'O-' w l, \
     'output.csv' u 54:19 title 'O_2-' w l, \
     'output.csv' u 54:20 title 'O_3-' w l, \
     'output.csv' u 54:21 title 'O_4-' w l 


set title "Negative Nitrogen-Oxide Ion Species"
set key left top
plot 'output.csv' u 54:22 title 'NO-' w l, \
     'output.csv' u 54:23 title 'N_2O-' w l, \
     'output.csv' u 54:24 title 'NO_2-' w l, \
     'output.csv' u 54:25 title 'NO_3-' w l 


set title "Negative Hydrogen Ion Species"
set key left top
plot 'output.csv' u 54:26 title 'H-' w l, \
     'output.csv' u 54:27 title 'OH-' w l, \



#NEUTRAL SPECIES ===================================

set title "Nitrogen Excited Species"
set key left top
plot 'output.csv' u 54:28 title 'N(2D)' w l, \
     'output.csv' u 54:29 title 'N_2(A3Sigma)' w l, \
     'output.csv' u 54:30 title 'N_2(B3Pi)' w l


set title "Atomic and Excited/Reactive Oxygen Species"
set key left top
plot 'output.csv' u 54:31 title 'O(1D)' w l, \
     'output.csv' u 54:32 title 'H' w l, \
     'output.csv' u 54:33 title 'N' w l, \
     'output.csv' u 54:34 title 'O' w l, \
     'output.csv' u 54:35 title 'O_2(a1Delta)' w l, \
     'output.csv' u 54:36 title 'O_3' w l



set title "Nitrogen-Oxide Species"
set key left top
plot 'output.csv' u 54:37 title 'NO' w l, \
     'output.csv' u 54:38 title 'N_2O' w l, \
     'output.csv' u 54:39 title 'NO_2' w l, \
     'output.csv' u 54:40 title 'NO_3' w l, \
     'output.csv' u 54:41 title 'N_2O_3' w l, \
     'output.csv' u 54:42 title 'N_2O_4' w l, \
     'output.csv' u 54:43 title 'N_2O_5' w l


set title "Hydrogen containing neutral Species"
set key left top
plot 'output.csv' u 54:44 title 'H_2' w l, \
     'output.csv' u 54:45 title 'OH' w l, \
     'output.csv' u 54:46 title 'HO_2' w l, \
     'output.csv' u 54:47 title 'H_2O_2' w l, \
     'output.csv' u 54:48 title 'HNO' w l, \
     'output.csv' u 54:49 title 'HNO_2' w l, \
     'output.csv' u 54:50 title 'HNO_3' w l


set yrange [1E19:1E26]
set title "N_2, O_2 and H_2O"
plot 'output.csv' u 54:51 title 'N_2' w l, \
     'output.csv' u 54:52 title 'O_2' w l, \
     'output.csv' u 54:53 title 'H_2O' w l
