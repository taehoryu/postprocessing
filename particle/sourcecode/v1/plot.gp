set terminal postscript eps enhanced color font 'H,25'



reset
set logscale y
set format y "10^{%T}"
set xlabel "t [s]" font "H,30"
set ylabel "height [km]" font "H,30"
set output 'tracer_particle_average.ps'
set yrange[10:10**5]
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

set style fill transparent solid 0.5
set label "radiation+pert. (512 cells)" at graph 0.03,0.30
set label "N_{vol}    =100" at graph 0.03,0.23
set label "A_{vol}/y_{vol}=10^{-3}" at graph 0.03,0.16
set label "{/Symbol s}_{vol}/y_{vol}=10^{-2}"at graph 0.03,0.09

!set key bottom right
plot \
"Timestamp_average" using 2:8:7 notitle "cell # 512" with filledcurves ls 10,\
"Timestamp_average" using 2:6 notitle "cell # 512" with lines ls 10,\
"Timestamp_average" using 2:7 notitle "cell # 512" with lines ls 20,\
"Timestamp_average" using 2:8 notitle "cell # 512" with lines ls 20,\



reset
set logscale y
set format y "10^{%T}"
set xlabel "t [s]" font "H,30"
set ylabel "height [km]" font "H,30"
set output 'tracer_particle_each.ps'
set yrange[10:10**5]
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

set style fill transparent solid 0.5
set label "radiation+pert. (512 cells)" at graph 0.03,0.30
set label "N_{vol}    =100" at graph 0.03,0.23
set label "A_{vol}/y_{vol}=10^{-3}" at graph 0.03,0.16
set label "{/Symbol s}_{vol}/y_{vol}=10^{-2}"at graph 0.03,0.09

!set key bottom right
plot for  \
[i=1:20] "Timestamp_parsed" using ($1==i?$4:1/0):($3/10**5) notitle "cell # 512" with linespoints ls i pi -10



