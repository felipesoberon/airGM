set terminal pdf
set output "densities.pdf"

set logscale y
set logscale x

set format x "%g"

set xlabel "Time (s)"
set ylabel "Density (m-3)"

set yrange [1E15:1E22]
set xrange [1E-8:1E3]
set grid
show grid

set datafile separator ","

set title "Species plotted in Kossyi et al., 1992"
set key right bottom
plot 'output.csv' u 54:34 title 'O' w l,      \
     'output.csv' u 54:36 title 'O_3' w l,    \
     'output.csv' u 54:37 title 'NO' w l,     \
     'output.csv' u 54:39 title 'NO_2' w l,   \
     'output.csv' u 54:43 title 'N_2O_5' w l, \
     'output.csv' u 54:38 title 'N_2O' w l,   \
     'output.csv' u 54:40 title 'NO_3' w l,   \
     'output.csv' u 54:33 title 'N' w l,      \
     'output.csv' u 54:42 title 'N_2O_4' w l
     




