######
# PHASE SPACE ERROR (MOVE FROM HERE)
#####

# reset
# set terminal pngcairo size 920,800 font 'Helvetica,15'
# set loadpath "../output"
# set output "../output/fig_orbital_angular_mentum.png"

# set xlabel "Cycle"
# set ylabel "Error"
# unset key
# set title "Orbital angular momentum error"

# plot 'phase_space_orbital_angular_momentum_error.dat' w d



######
# BASIN OF ATTRACTION
#####

reset
set terminal pngcairo size 920,800 font 'Helvetica,15'
set loadpath "../output"
set output "../output/fig_basin_of_attraction.png"

# set xrange[0.0:6.29]
set xrange[-3.1415:3.1415]
set yrange[0.0:3.0]
unset key
unset title

plot 'basin.dat' u 1:2:3 w image notitle

reset
set terminal pngcairo size 920,800 font 'Helvetica,15'
set loadpath "../output"
set output "../output/fig_basin_of_attraction_times.png"

# set xrange[0.0:6.29]
set xrange[-3.1415:3.1415]
set yrange[0.0:3.0]
unset key
unset title

plot 'basin.dat' u 1:2:4 w image notitle

reset
set terminal pngcairo size 920,800 font 'Helvetica,15'
set loadpath "../output"
set output "../output/fig_basin_of_attraction_rotation.png"
unset key
unset title

# set xrange[0.0:6.28318530718]
set xrange[-3.1415:3.1415]
set yrange[0.0:3.0]
unset key

plot 'basin.dat' u 1:2:5 w image notitle

######
# BASIN OF ATTRACTION - UNION
#####

reset
set terminal pngcairo size 920,800 font 'Helvetica,15'
set loadpath "../output"
set output "../output/fig_basin_of_attraction_union.png"

# set xrange[0.0:6.29]
set xrange[-3.1415:3.1415]
set yrange[0.0:3.0]
unset key
unset title

plot 'basin_union.dat' u 1:2:3 w image notitle

reset
set terminal pngcairo size 920,800 font 'Helvetica,15'
set loadpath "../output"
set output "../output/fig_basin_of_attraction_union_times.png"

# set xrange[0.0:6.29]
set xrange[-3.1415:3.1415]
set yrange[0.0:3.0]
unset key
unset title

plot 'basin_union.dat' u 1:2:4 w image notitle

reset
set terminal pngcairo size 920,800 font 'Helvetica,15'
set loadpath "../output"
set output "../output/fig_basin_of_attraction_union_rotation.png"
unset key
unset title

# set xrange[0.0:6.28318530718]
set xrange[-3.1415:3.1415]
set yrange[0.0:3.0]
unset key

plot 'basin_union.dat' u 1:2:5 w image notitle