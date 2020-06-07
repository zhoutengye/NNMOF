set log y
set key right bottom
set yrange [6.28319e-4:1]
# Compute growth rate assuming ∆U = 1
#2 π = stemp=`awk -v "r=1" -v "w=1" 'BEGIN {ak=2.*w*atan2(0, -1); print sqrt(r)*ak}'`
plot "stats" u 1:16 t "r = 1" w l, 6.28319e-04*exp(6.28319*x) t "theory r=1" 


