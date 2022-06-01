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

# reset
# set terminal pngcairo size 920,800 font 'Helvetica,15'
# set loadpath "../output"
# set output "../output/fig_basin_of_attraction.png"

# # set xrange[0.0:6.29]
# set xrange[-3.1415:3.1415]
# set yrange[0.0:3.0]
# unset key
# set title "Phase space {/Symbol q}-~{/Symbol q}{1.1.}"

# plot 'basin.dat' u 1:2:3 w image

# reset
# set terminal pngcairo size 920,800 font 'Helvetica,15'
# set loadpath "../output"
# set output "../output/fig_basin_of_attraction_times.png"

# # set xrange[0.0:6.29]
# set xrange[-3.1415:3.1415]
# set yrange[0.0:3.0]
# unset key
# set title "Phase space {/Symbol q}-~{/Symbol q}{1.1.}"

# plot 'basin.dat' u 1:2:4 w image

reset
set terminal pngcairo size 920,800 font 'Helvetica,15'
set loadpath "../output"
set output "../output/fig_basin_of_attraction_union.png"

# set xrange[0.0:6.29]
set xrange[-3.1415:3.1415]
set yrange[0.0:3.0]
unset key
set title "Phase space {/Symbol q}-~{/Symbol q}{1.1.}"

plot 'basin_union.dat' u 1:2:3 w image

reset
set terminal pngcairo size 920,800 font 'Helvetica,15'
set loadpath "../output"
set output "../output/fig_basin_of_attraction_union_times.png"

# set xrange[0.0:6.29]
set xrange[-3.1415:3.1415]
set yrange[0.0:3.0]
unset key
set title "Phase space {/Symbol q}-~{/Symbol q}{1.1.}"

plot 'basin_union.dat' u 1:2:4 w image