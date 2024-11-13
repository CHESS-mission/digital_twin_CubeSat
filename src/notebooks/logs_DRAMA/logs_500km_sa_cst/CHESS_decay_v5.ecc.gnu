reset
set term png
set output "CHESS_decay_v5.ecc.png"
set xdata time
set timefmt "%Y-%m-%d %H:%M:%S"
set format x "%Y/%m"
set xtics 7884000
set title "DRAMA\nOSCAR - Orbital Spacecraft Active Removal\nEccentricity vs. Time"
set xlabel "Date"
set ylabel "Eccentricity\nSingly averaged (over M)"
set xrange [*:*]
set yrange [*:*]
set format y "%7.5f"
set key below
plot \
"CHESS_decay_v5.oev" u 1:5 w l lt 01 lw 01 notitle
