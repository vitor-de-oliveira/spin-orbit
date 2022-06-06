reset
set terminal pngcairo size 920,800 font 'Helvetica,15'
set loadpath "../output"
set output "../output/fig_time_series.png"

unset key
unset title

plot for [i=1:10] 'orbit_'.i.'.dat' u 2 w l lw 2

reset
set terminal pngcairo size 920,800 font 'Helvetica,15'
set loadpath "../output"
set output "../output/fig_time_series_zoom.png"

set yrange[0:10]
unset key
unset title

plot for [i=1:10] 'orbit_'.i.'.dat' u 2 w l lw 2

reset
set terminal pngcairo size 920,800 font 'Helvetica,15'
set loadpath "../output"
set output "../output/fig_time_series_log.png"

set log y
unset key
unset title

plot for [i=1:10] 'orbit_'.i.'.dat' u 2 w l lw 2
