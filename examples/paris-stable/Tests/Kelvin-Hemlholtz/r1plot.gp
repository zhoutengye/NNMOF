set log y
set format y "% g"
set key right bottom
set xlabel "t"
set ylabel "maximum |v|"
plot 'stats' u 1:16 t "num", 1e-2*exp(3.14*x) t 'Theory s=2π'

