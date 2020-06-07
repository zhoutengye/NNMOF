set multiplot layout 2,3
set log y
set format y "% g"
plot "out/ab_coef_file_.txt" u 1:4 w l t "norm of (a,b)"
set key right bottom
plot 'stats' u 1:16 t "|v|"

