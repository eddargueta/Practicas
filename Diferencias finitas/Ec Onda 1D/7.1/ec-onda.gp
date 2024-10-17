# Cargar este script utilizando gnuplot con: load 'ec-onda.gp'.

set xrange [0:1]
set yrange [-2:2]

set xlabel 'x [m]'
set ylabel 'u(x,t)'

set xtics 0.25

dt=0.0001

do for [it=0:200] {
    set title sprintf("t = %f (s)", it*dt)
    plot 'solucion-1.dat' index it u 2:3 w l t 'Solucion', 
    # plot 'solucion-2.dat' index it u 2:3 w l t 'Solucion', 
    # plot 'solucion-3.dat' index it u 2:3 w l t 'Solucion',
    # plot 'solucion-4.dat' index it u 2:3 w l t 'Solucion' 
    pause -1 # Presionar enter para ver el movimiento de la onda.
    # pause 0.05 # Para que la animacion sea en automatico. 
    replot
}







