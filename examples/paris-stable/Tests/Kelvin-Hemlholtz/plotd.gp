set multiplot layout 3,2
set log y
set format y "% g"
set key right bottom
plot "out/ab_coef_file_.txt" u 1:4 w l t "modulus of Fourier ampl. for wnr 1"
plot 'stats' u 1:16 t "max of |v|"
#
plot "out/ab_coef_file_.txt" u 1:6 w l t "modulus of Fourier amplitude of wnr 2"
plot "out/ab_coef_file_.txt" u 1:7 w l t "modulus of Fourier amplitude of wnr 3"
#
plot "out/ab_coef_file_.txt" u 1:8 w l t "modulus of Fourier amplitude of wnr 4"
plot "out/ab_coef_file_.txt" u 1:9 w l t "modulus of Fourier amplitude of wnr 5"
