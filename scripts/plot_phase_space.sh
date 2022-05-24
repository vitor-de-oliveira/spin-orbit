# reset
# set terminal pngcairo size 920,800 font 'Helvetica,15'
# set loadpath "../output"
# set output "../output/fig_orbital_angular_mentum.png"

# set xlabel "Cycle"
# set ylabel "Error"
# unset key
# set title "Orbital angular momentum error"

# plot 'phase_space_orbital_angular_momentum_error.dat' w d

# reset
# set terminal pngcairo size 920,800 font 'Helvetica,15'
# set loadpath "../output"
# set output "../output/fig_phase_space.png"

# set xrange[0.0:6.29]
# set yrange[0.0:3.0]
# unset key
# set title "Phase space {/Symbol q}-~{/Symbol q}{1.1.}"

# plot 'phase_space.dat' w d, \
#      # 'phase_space_initial_conditions.dat' w p, \

# reset
# set terminal pngcairo size 920,800 font 'Helvetica,15'
# set loadpath "../output"
# set output "../output/fig_phase_space_with_dissipative_orbit.png"

# set xrange[0.0:6.29]
# set yrange[0.0:3.0]
# unset key
# set title "Phase space {/Symbol q}-~{/Symbol q}{1.1.}"

# plot 'phase_space.dat' w d, \
#      'orbit.dat' w p pt 7 ps 1.5

# reset
# set terminal pngcairo size 920,800 font 'Helvetica,15'
# set loadpath "../output"
# set output "../output/fig_phase_space_with_dissipative_orbit_zoom.png"

# set xrange[2.45:3.8]
# set yrange[0.1:1.2]
# unset key
# set title "Phase space {/Symbol q}-~{/Symbol q}{1.1.}"

# plot 'phase_space.dat' w d, \
#      'orbit.dat' w p pt 7 ps 1.5

# reset
# set terminal pngcairo size 920,800 font 'Helvetica,15'
# set loadpath "../output"
# set output "../output/fig_basin_of_attraction.png"

# set xrange[0.0:6.29]
# set yrange[0.0:3.0]
# unset key
# set title "Phase space {/Symbol q}-~{/Symbol q}{1.1.}"

# plot 'phase_space.dat' w d, \
#      'basin.dat' u 1:2:3 w image, \

reset
set terminal pngcairo size 920,800 font 'Helvetica,15'
set loadpath "../output"
set output "../output/fig_basin_of_attraction.png"

# set xrange[0.0:6.28318530718]
set xrange[-3.1415:3.1415]
set yrange[0.0:3.0]
unset key
set title "Phase space {/Symbol q}-~{/Symbol q}{1.1.}"

# plot 'phase_space.dat' w d, \
#      'basin.dat' u 1:2:3 w image, \

plot 'basin.dat' u 1:2:3 w image, \

# reset
# set terminal pngcairo size 920,800 font 'Helvetica,15'
# set loadpath "../output"
# set output "../output/fig_basin_of_attraction_no_opt.png"

# # set xrange[0.0:6.28318530718]
# set xrange[-3.1415:3.1415]
# set yrange[0.0:3.0]
# unset key
# set title "Phase space {/Symbol q}-~{/Symbol q}{1.1.}"

# # plot 'phase_space.dat' w d, \
# #      'basin.dat' u 1:2:3 w image, \

# plot 'basin_no_opt.dat' u 1:2:3 w image, \

# reset
# set terminal pngcairo size 920,800 font 'Helvetica,15'
# set loadpath "../output"
# set output "../output/fig_basin_of_attraction_union.png"

# # set xrange[0.0:6.28318530718]
# set xrange[-3.1415:3.1415]
# set yrange[0.0:3.0]
# unset key
# set title "Phase space {/Symbol q}-~{/Symbol q}{1.1.}"

# # plot 'phase_space.dat' w d, \
# #      'basin.dat' u 1:2:3 w image, \

# plot 'basin_union.dat' u 1:2:3 w image, \

# reset
# set terminal pngcairo size 920,800 font 'Helvetica,15'
# set loadpath "../output"
# set output "../output/fig_basin_of_attraction_no_omp.png"

# # set xrange[0.0:6.28318530718]
# set xrange[-3.1415:3.1415]
# set yrange[0.0:3.0]
# unset key
# set title "Phase space {/Symbol q}-~{/Symbol q}{1.1.}"

# # plot 'phase_space.dat' w d, \
# #      'basin.dat' u 1:2:3 w image, \

# plot 'basin_no_omp.dat' u 1:2:3 w image, \

# reset
# set terminal pngcairo size 920,800 font 'Helvetica,15'
# set loadpath "../output"
# set output "../output/fig_basin_of_attraction_no_grid.png"

# # set xrange[0.0:6.28318530718]
# set xrange[-3.1415:3.1415]
# set yrange[0.0:3.0]
# unset key
# set title "Phase space {/Symbol q}-~{/Symbol q}{1.1.}"

# # plot 'phase_space.dat' w d, \
# #      'basin.dat' u 1:2:3 w image, \

# plot 'basin_no_grid.dat' w p pt 7 ps 1

# reset
# set terminal pngcairo size 920,800 font 'Helvetica,15'
# set loadpath "../output"
# set output "../output/fig_basin_of_attraction_no_grid_no_omp.png"

# # set xrange[0.0:6.28318530718]
# set xrange[-3.1415:3.1415]
# set yrange[0.0:3.0]
# unset key
# set title "Phase space {/Symbol q}-~{/Symbol q}{1.1.}"

# # plot 'phase_space.dat' w d, \
# #      'basin.dat' u 1:2:3 w image, \

# plot 'basin_no_grid_no_omp.dat' w p pt 7 ps 1

# reset
# set terminal pngcairo size 920,800 font 'Helvetica,15'
# set loadpath "../output"
# set output "../output/fig_basin_of_attraction_no_omp.png"

# # set xrange[0.0:6.28318530718]
# set xrange[-3.1415:3.1415]
# set yrange[0.0:3.0]
# unset key
# set title "Phase space {/Symbol q}-~{/Symbol q}{1.1.}"

# # plot 'phase_space.dat' w d, \
# #      'basin.dat' u 1:2:3 w image, \

# plot 'basin_no_omp.dat' u 1:2:3 w image, \

# reset
# set terminal pngcairo size 920,800 font 'Helvetica,15'
# set loadpath "../output"
# set output "../output/fig_basin_of_attraction_no_omp.png"

# # set xrange[0.0:6.28318530718]
# set xrange[-3.1415:3.1415]
# set yrange[0.0:3.0]
# unset key
# set title "Phase space {/Symbol q}-~{/Symbol q}{1.1.}"

# plot 'basin_no_omp.dat' u 1:2:3 w image, \

# reset
# set terminal pngcairo size 920,800 font 'Helvetica,15'
# set loadpath "../output"
# set output "../output/fig_basin_of_attraction_times.png"

# set xrange[0.0:6.29]
# set yrange[0.0:3.0]
# unset key
# set title "Phase space {/Symbol q}-~{/Symbol q}{1.1.}"

# plot 'basin.dat' u 1:2:4 w image, \