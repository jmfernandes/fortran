set terminal epslatex color
set output "eg2.tex"
set format y "$%g$"
set format x "$%g$"
set title 'double pendulum - low energy'
set xlabel '$\theta_2$'
set ylabel '$p_{\theta_2}$'

set grid lt 1 lw 1 lc rgb "#becc36"

set xrange [-0.025:0.025]
#set yrange [-0.025:0.025]

set xtic -0.025,0.01,0.025
#set ytic -0.025,0.01,0.025

set key top right

#set ytics ("0" -2e11, "0.5" 0, "1" 2e11) nomirror
#set xtics ('$\pi$' -2e11,\
#'$x$' 0,\
#'$hey$' 2e11)
plot "angle_mom2.data" using 1:2 with lines linewidth 1 linetype 1 linecolor 3 title "pendulum 2"

set output
set terminal pop