set terminal postscript eps enhanced color font 'H,20'


reset

set xlabel "x [km]" font "H,20"
set ylabel "y [km]" font "H,20"
set output 'tracer_trajectory_1mbar.ps'
set style line 10 linewidth 3.000 dashtype 1 pointtype 7 pointsize 1.0 pointinterval 0 lc 7
set style line 20 linewidth 1.000 dashtype 1 pointtype 7 pointsize 1.0 pointinterval 0 lc 7
set style line 12 linewidth 1.000 dashtype 1 pointtype 4 pointsize 1.0 pointinterval 0 lc rgb"blue"
set style line 40 linewidth 2.000 dashtype 1 pointtype 1 pointsize 2.0 pointinterval 0 lc "#008000"


set style line 21 linewidth 10.000 dashtype 1 pointtype 6 pointsize 3.0 pointinterval 0 lc rgb"red"
set style line 22 linewidth 6.000 dashtype 1 pointtype 6 pointsize 2.0 pointinterval 0 lc rgb"red"
set style line 23 linewidth 4.000 dashtype 1 pointtype 6 pointsize 1.5 pointinterval 0 lc rgb"red"

set style line 31 linewidth 10.000 dashtype 1 pointtype 7 pointsize 3.0 pointinterval 0 lc rgb"blue"
set style line 32 linewidth 6.000 dashtype 1 pointtype 7 pointsize 2.0 pointinterval 0 lc rgb"blue"
set style line 33 linewidth 4.000 dashtype 1 pointtype 7 pointsize 1.5 pointinterval 0 lc rgb"blue"

set label "hydro+pert. (2048 cells)" at graph 0.03,0.30
set label "N_{vol}    =10" at graph 0.03,0.23
set label "A_{vol}/y_{vol}=10^{-3}" at graph 0.03,0.16
set label "{/Symbol s}_{vol}/y_{vol}=10^{-2}"at graph 0.03,0.09
set yrange[0:1.5*10**4]
set multiplot
set xrange[0:2.25*10**4]
plot for  \
[i=1:500] 'P001mbar/trajectory_'.i.'.txt' using ($2/10**5):($3/10**5) notitle "cell # 512" with points ls i
plot for  \
[i=1:500] 'P010mbar/trajectory_'.i.'.txt' using ($2/10**5):($3/10**5) notitle "cell # 512" with points ls i

plot for  \
[i=1:500] 'P100mbar/trajectory_'.i.'.txt' using ($2/10**5):($3/10**5) notitle "cell # 512" with points ls i
plot for  \
[i=1:500] 'P0001bar/trajectory_'.i.'.txt' using ($2/10**5):($3/10**5) notitle "cell # 512" with points ls i
plot for  \
[i=1:500] 'P0010bar/trajectory_'.i.'.txt' using ($2/10**5):($3/10**5) notitle "cell # 512" with points ls i


unset multiplot




