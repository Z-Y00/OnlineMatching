set output "6.png"
set terminal png truecolor
set autoscale
#set xrange [0:20]
#set yrange [0:20]

set style increment default
set style arrow 1 head back filled linecolor rgb "dark-violet"  linewidth 2.000 dashtype 1 size screen  0.025,30.000,45.000
set style arrow 2 head back nofilled linecolor rgb "skyblue"  linewidth 2.000 dashtype 1 size screen  0.030,15.000,90.000
set style arrow 3 head back filled linecolor rgb "dark-violet"  linewidth 2.000 dashtype 1 size screen  0.030,15.000,45.000
myencoding = "utf8"

plot 'servers.dat', 'out.dat'  with vectors arrowstyle 3
