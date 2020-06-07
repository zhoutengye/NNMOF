set xrange [0:0.025]
set xlabel "time"
set ylabel "E_k"
set yrange [*:*]
set title "Oscillating bubble simulation with VOF"
plot "reference_Ek.txt" t "D/∆x=38 reference simulation" w l, "test_Ek.txt" t "D/∆x=19 test simulation" w p, 0 notitle
# Remove because of OS 10.10 problem
#set term pdf
#set out "tmp.pdf"
#replot


