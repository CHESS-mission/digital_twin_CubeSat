reset
set term png
set output "CHESS_decay_v9.alt.png"
set xdata time
set timefmt "%Y-%m-%d %H:%M:%S"
set format x "%Y"
set title "DRAMA\nOSCAR - Orbital Spacecraft Active Removal\nAltitude vs. Time"
set xlabel "Date"
set xrange [*:*]
set yrange [*:*]
set ylabel  "Altitude [km]\nSingly averaged (over M)"
set key below
plot \
"CHESS_decay_v9.oev" u 1:($4*(1.0 - $5) - 6378.137) w l lt 01 lw 01 title "Perigee altitude", \
"CHESS_decay_v9.oev" u 1:($4*(1.0 + $5) - 6378.137) w l axes x1y1 lt 02 lw 01 title "Apogee altitude"
