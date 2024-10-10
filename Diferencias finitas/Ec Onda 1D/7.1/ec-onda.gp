set xrange [0:1]
set yrange [-1:1]

set xlabel 'x [m]'
set ylabel 'u(x,t)'

dt=0.0001

do for [it=0:200] {
    set title sprintf("t = %f (s)", it*dt)
    plot 'solucion-1d-1.dat' index it u 2:3 w l t 'Solucion'
    
    pause -1
}







