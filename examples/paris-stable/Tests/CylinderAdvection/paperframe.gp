#set term png
#set output 'bubble.png'

#set xrange [0:1]
#set yrange [0:1]

set size square
unset ytics
unset xtics

#set xlabel " x/h "
#set ylabel " y/h "
set term pngcairo
set out "tmp.png"
unset surface
set contour
set cntrparam linear
set view map
unset clabel
set cntrparam levels discrete 0.5

#set title "ouptut nr  ".frame
#set grid

splot  "toplot1.tmp" lt -1 w l   notitle  , "toplot2.tmp" lt -1 w l   notitle  , "toplot3.tmp" lt -1 w l   notitle   , "toplot4.tmp" lt -1 w l   notitle 
exit


