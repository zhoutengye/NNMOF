set xrange [0.075:*]
set log y
#plot "stats" u 1:13 t "K energy phase C=1 rho2 liq", "stats" u 1:14 t "K energy phase C=0 rho1 gas"
plot  "stats" u 1:14 t "K energy phase C=0 rho1 gas" w l

