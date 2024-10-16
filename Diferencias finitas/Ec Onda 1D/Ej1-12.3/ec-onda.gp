set xrange [0:1]
set yrange [-2:2]

set xlabel 'x [m]'
set ylabel 'u(x,t)'

set xtics 0.25

dt=0.001

do for [it=0:400] {
    set title sprintf("t = %f (s)", it*dt)
    plot 'solucion-1d-1.dat' index it u 2:3 w l t 'alpha = 0.25', 'solucion-1d-2.dat' index it u 2:3 w l t 'alpha = 0.50', 'solucion-1d-3.dat' index it u 2:3 w l t 'alpha = 0.75', 'solucion-1d-4.dat' index it u 2:3 w l t 'alpha = 1.00', 'solucion-1d-5.dat' index it u 2:3 w l t 'alpha = 1.25'	
    pause -1
}






