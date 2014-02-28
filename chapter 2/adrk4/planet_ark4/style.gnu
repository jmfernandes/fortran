set terminal aqua enhanced font "arial,15" size 900, 600
set output 'plot.png'
set grid
set key top right
set border lw 1 lc rgb "black"



# change text colors of tics
set xtics textcolor rgb "black" offset -1
set ytics textcolor rgb "black" offset 3
 
# change text colors of labels
set xlabel "X" textcolor rgb "black"
set ylabel "Y" textcolor rgb "black"
 
# change a text color of key
set key textcolor rgb "black"

#set grid color
set grid ytics lt 0 lw 1 lc rgb "#F5B800"


#######################################################################

set view 45,45
set contour base

set title 'Figure'
set xlabel 'Distance (meters)'
set ylabel 'Distance (meters)'
set zlabel 'Distance (meters)' offset -3 rotate left

#for auto color
splot 'data.data' using 2:3:4 with lines title 'data 1',\
'data.data' using 2:3:(0) with lines title 'data 2'
#'data.dat' using 1:3 with lines title 'data 2',\
#'data.dat' using 1:4 with lines title 'data 3'

#to specify the color exaclty
#plot 'data.dat' using 1:2 with linespoints ls 0 lc 3 title 'data 1',\
#'data.dat' using 1:3 with linespoints ls 0 lc 4 title 'data 2',\
#'data.dat' using 1:4 with linespoints ls 0 lc 5 title 'data 3'




#single plot
#plot 'data.dat' title 'data 1' with linespoints ls 0 lc 3

#for multiple plots use this
#plot 'data.dat' using 1:2:(column(-2)) with lines lc variable title 'data 1'

# ls
# 0= no points
# 1= crosses
# 2= x's
# 3= asterisk
# 4= squares
# 5= diamonds
# 6= triangles

#lc color
# 0=grey
# 1=red
# 2=bright green
# 3=blue
# 4=pink
# 5=teal
# 6=brown
# 7=orange
# 8=brownish green
# 9=darker brownish green