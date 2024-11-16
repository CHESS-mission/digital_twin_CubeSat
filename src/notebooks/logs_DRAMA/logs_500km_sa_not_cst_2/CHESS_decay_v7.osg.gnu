reset
set term png
set output "CHESS_decay_v7.osg.png"
set xdata time
set timefmt "%d/%m/%Y"
set format x "%Y"
set title "DRAMA\nOSCAR - Orbital Spacecraft Active Removal\nSolar activity vs. Time"
set xlabel "Date"
set ylabel "Solar activity [sfu]"
set xrange [*:*]
set yrange [*:*]
set format y "%5.1f"
set key below
plot \
"CHESS_decay_v7.osg" u 1:2 w p lt 01 lw 01 title "Daily F10.7", \
"CHESS_decay_v7.osg" u 1:3 w l lt 02 lw 01 title "Avg. F10.7 (Last 81 days)"
