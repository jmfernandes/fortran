set terminal epslatex color
set output "eg1.tex"
set format y "$%g$"
set format x "$%g$"
set title 'single pendulum'
set xlabel '$\theta$'
set ylabel '$p_{\theta}$'

set grid lt 1 lw 1 lc rgb "#becc36"

#set xrange [-2e11:2e11]
#set yrange [-2e11:2e11]

#set xtic -2e11,1e11,2e11
#set ytic -2e11,0.5e11,2e11

set key top right

#set ytics ("0" -2e11, "0.5" 0, "1" 2e11) nomirror
#set xtics ('$\pi$' -2e11,\
#'$x$' 0,\
#'$hey$' 2e11)
plot "aapend_ang_mo.data" using 1:2 with lines linewidth 1 linetype 1 linecolor 3 title "plot"

set output
set terminal pop