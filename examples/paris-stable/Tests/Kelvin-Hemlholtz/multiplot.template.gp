set log y
set format y "% g"
set key right bottom
set yrange [*:*]
set title 'TITLE'
plot "SRHOG" u 1:16 t "r = RHOG",  "SRHOG" u 1:16 t "r = RHOG",  "SRHOG" u 1:16 t "r = RHOG",  "SRHOG" u 1:16 t "r = RHOG",  "SRHOG" u 1:16 t "r = RHOG", ATEMP*exp(STEMP*x) t "theory r=DRHOG", ATEMP*exp(STEMP*x) t "theory r=DRHOG"
