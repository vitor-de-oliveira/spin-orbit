reset
set terminal pngcairo size 920,800 font 'Helvetica,15'
set loadpath "../output"
set output "../output/fig_histogram_times_union.png"

n=50 #number of intervals
max=450.
min=0.
width=(max-min)/n 
hist(x,width)=width*floor(x/width)+width/2.0
set xrange [min:max]
set yrange [0:]
# set offset graph 0.05,0.05,0.05,0.0
set xtics min,(max-min)/5,max
set boxwidth width*0.9
set style fill solid
set tics out nomirror
set xlabel "x"
set ylabel "Frequency"

# set table '../output/hist.dat'
# plot "basin.dat" u (hist($4,width)):(1.0) smooth freq w boxes notitle
# unset table

# gauss(x)=a/(sigma*sqrt(2.*pi))*exp(-(x-mu)**2./(2.*sigma**2))
# a=1000; sigma = 100; mu = -100
# fit gauss(x) 'hist.temp' u 1:2 via a, sigma, mu

# f(x) = b + n*exp(-x/u)
# fit[70:][] log(f(x)) "hist.temp" using 1:(log($2)) via b,n,u
# fit[70:100][0:60] f(x) "hist.temp" using 1:2 via b,n,u

plot "basin_union.dat" u (hist($4,width)):(1.0) smooth freq w boxes notitle
# , f(x) w lines ls 2 lw 2