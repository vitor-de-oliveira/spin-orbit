reset
set terminal pngcairo size 920,800 font 'Helvetica,15'
set loadpath "../output/basin_of_attraction/"
set output "../output/tests/fig_multiple_basins.png"

set xrange [-3.15:3.15]
set yrange [0:3]

set key opaque

unset colorbox

res = 600

plot 'multiple_basin_determined_gamma_0.264_e_0.100_system_linear_K_0.01000_res_'.res.'_n_1000_basin_eps_0.100.dat' u 1:2:($3==1||$3==3?1:1/0) w image notitle, \
     'multiple_basin_determined_gamma_0.264_e_0.100_system_linear_K_0.01000_res_'.res.'_n_1000_basin_eps_0.100.dat' u 1:2:($3==2?2:1/0) w image notitle, \
     'multiple_basin_determined_ref_gamma_0.264_e_0.100_system_linear_K_0.01000_res_'.res.'_n_1000_basin_eps_0.100.dat' u 1:2 w p pt 7 ps 2 lc rgb "red" title "spin-orbit resonances"