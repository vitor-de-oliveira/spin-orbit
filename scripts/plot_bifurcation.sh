reset
set terminal pngcairo size 920,800 font 'Helvetica,20'
set loadpath "../output/periodic_orbit"

set key box opaque left
set xlabel "e"
set ylabel "|{/Symbol l}|"

# set output "../output/meeting_10_21/bifurcation_resonance_1_1.png"

# set title "Resonance 1/1"

# plot 'periodic_orbit_track_gamma_0.264_system_rigid_resonance_1_1.dat' u 1:2 w l lw 2 lc rgb "black" title "K = 0.000", \
#      'periodic_orbit_track_gamma_0.264_system_rigid_resonance_1_1.dat' u 1:3 w l lw 2 lc rgb "black" notitle, \
#      'periodic_orbit_track_gamma_0.264_system_linear_K_0.00100_resonance_1_1.dat' u 1:2 w l lw 2 lc rgb "green" title "K = 0.001", \
#      'periodic_orbit_track_gamma_0.264_system_linear_K_0.00100_resonance_1_1.dat' u 1:3 w l lw 2 lc rgb "green" notitle, \
#      'periodic_orbit_track_gamma_0.264_system_linear_K_0.01000_resonance_1_1.dat' u 1:2 w l lw 2 lc rgb "red" title "K = 0.010", \
#      'periodic_orbit_track_gamma_0.264_system_linear_K_0.01000_resonance_1_1.dat' u 1:3 w l lw 2 lc rgb "red" notitle, \
#      'periodic_orbit_track_gamma_0.264_system_linear_K_0.05000_resonance_1_1.dat' u 1:2 w l lw 2 lc rgb "purple" title "K = 0.050", \
#      'periodic_orbit_track_gamma_0.264_system_linear_K_0.05000_resonance_1_1.dat' u 1:3 w l lw 2 lc rgb "purple" notitle, \
#      'periodic_orbit_track_gamma_0.264_system_linear_K_0.10000_resonance_1_1.dat' u 1:2 w l lw 2 lc rgb "blue" title "K = 0.100", \
#      'periodic_orbit_track_gamma_0.264_system_linear_K_0.10000_resonance_1_1.dat' u 1:3 w l lw 2 lc rgb "blue" notitle

# set output "../output/meeting_10_21/bifurcation_resonance_2_1.png"

# set title "Resonance 2/1"

# plot 'periodic_orbit_track_gamma_0.264_system_rigid_resonance_2_1.dat' u 1:2 w l lw 2 lc rgb "black" title "K = 0.000", \
#      'periodic_orbit_track_gamma_0.264_system_rigid_resonance_2_1.dat' u 1:3 w l lw 2 lc rgb "black" notitle, \
#      'periodic_orbit_track_gamma_0.264_system_linear_K_0.00100_resonance_2_1.dat' u 1:2 w l lw 2 lc rgb "green" title "K = 0.001", \
#      'periodic_orbit_track_gamma_0.264_system_linear_K_0.00100_resonance_2_1.dat' u 1:3 w l lw 2 lc rgb "green" notitle, \
#      'periodic_orbit_track_gamma_0.264_system_linear_K_0.01000_resonance_2_1.dat' u 1:2 w l lw 2 lc rgb "red" title "K = 0.010", \
#      'periodic_orbit_track_gamma_0.264_system_linear_K_0.01000_resonance_2_1.dat' u 1:3 w l lw 2 lc rgb "red" notitle

# set output "../output/meeting_10_21/bifurcation_resonance_1_2.png"

# set key box opaque bottom left

# set title "Resonance 1/2"

# plot 'periodic_orbit_track_gamma_0.264_system_rigid_resonance_1_2.dat' u 1:2 w l lw 2 lc rgb "black" title "K = 0.000", \
#      'periodic_orbit_track_gamma_0.264_system_rigid_resonance_1_2.dat' u 1:3 w l lw 2 lc rgb "black" notitle, \
#      'periodic_orbit_track_gamma_0.264_system_linear_K_0.00100_resonance_1_2.dat' u 1:2 w l lw 2 lc rgb "green" title "K = 0.001", \
#      'periodic_orbit_track_gamma_0.264_system_linear_K_0.00100_resonance_1_2.dat' u 1:3 w l lw 2 lc rgb "green" notitle, \
#      'periodic_orbit_track_gamma_0.264_system_linear_K_0.01000_resonance_1_2.dat' u 1:2 w l lw 2 lc rgb "red" title "K = 0.010", \
#      'periodic_orbit_track_gamma_0.264_system_linear_K_0.01000_resonance_1_2.dat' u 1:3 w l lw 2 lc rgb "red" notitle

set output "../output/meeting_10_21/bifurcation_resonance_3_2.png"

set key box opaque bottom left

set title "Resonance 3/2"

plot 'periodic_orbit_track_gamma_0.264_system_rigid_resonance_3_2.dat' u 1:2 w l lw 2 lc rgb "black" title "K = 0.000", \
     'periodic_orbit_track_gamma_0.264_system_rigid_resonance_3_2.dat' u 1:3 w l lw 2 lc rgb "black" notitle, \
     'periodic_orbit_track_gamma_0.264_system_linear_K_0.00100_resonance_3_2.dat' u 1:2 w l lw 2 lc rgb "green" title "K = 0.001", \
     'periodic_orbit_track_gamma_0.264_system_linear_K_0.00100_resonance_3_2.dat' u 1:3 w l lw 2 lc rgb "green" notitle, \
     'periodic_orbit_track_gamma_0.264_system_linear_K_0.01000_resonance_3_2.dat' u 1:2 w l lw 2 lc rgb "red" title "K = 0.010", \
     'periodic_orbit_track_gamma_0.264_system_linear_K_0.01000_resonance_3_2.dat' u 1:3 w l lw 2 lc rgb "red" notitle