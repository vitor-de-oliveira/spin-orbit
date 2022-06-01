reset
set terminal pngcairo size 920,800 font 'Helvetica,15'
set loadpath "../output"
set output "../output/fig_histogram_times_250.png"

n=200 #number of intervals
max=450. #max value
min=0. #min value
width=(max-min)/n #interval width
#function used to map a value to the intervals
hist(x,width)=width*floor(x/width)+width/2.0
set xrange [min:max]
set yrange [0:]
#to put an empty boundary around the
#data inside an autoscaled graph.
# set offset graph 0.05,0.05,0.05,0.0
set xtics min,(max-min)/5,max
set boxwidth width*0.9
set style fill solid #fillstyle
set tics out nomirror
set xlabel "x"
set ylabel "Frequency"

gauss(x)=a/(sigma*sqrt(2.*pi))*exp(-(x-mu)**2./(2.*sigma**2))

fit gauss(x) "basin_250.dat" u (hist($4,width)):(1.0) via a, sigma, mu

#count and plot
plot "basin_250.dat" u (hist($4,width)):(1.0) smooth freq w boxes lc rgb "green" notitle