import os
import argparse
import numpy as np
import scipy.stats as st
import statsmodels as sm
from matplotlib import cm
import matplotlib.pyplot as plt
from fractions import Fraction as frac
from matplotlib.ticker import FuncFormatter
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import FormatStrFormatter

def __get_arg() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument('--input',
                        action = 'store',
                        dest = 'input_path')
    parser.add_argument('--output',
                        action = 'store',
                        dest = 'output_path')
    parser.add_argument('--parameter',
                        action = 'store',
                        dest = 'parameter_value') 
    return parser

def get_input() -> str:
    parser = __get_arg()
    args = parser.parse_args()
    input_path = args.input_path
    return input_path

def get_output() -> str:
    parser = __get_arg()
    args = parser.parse_args()
    output_path = args.output_path
    return output_path

def get_parameter() -> str:
    parser = __get_arg()
    args = parser.parse_args()
    par = args.parameter_value
    return par

def data_parse(flname):
    data = []
    with open(flname, 'r') as inpt:
        for line in inpt:
            if line != '' and line != '\n':
                data.append(float(line.split()[3]))
    return data

def plot_histogram(data, Nbins, output_name, par_value):
    plt.rc('text', usetex = True)  # LaTeX text
    plt.rc('font', family = 'serif', size = 70)
    fig, ax = plt.subplots(figsize = (30, 28))

    # Set X axes of histogram #
    # Xmin = min(data)
    # Xmax = max(data)
    # XMt = MultipleLocator(abs(Xmax - Xmin) / 5.00)  # X major ticks
    # Xmt = MultipleLocator(abs(Xmax - Xmin) / 20.0)  # X minor ticks
    # X_format = FormatStrFormatter('%0.1f')          # decimal labels
    # ax.set_xlim([Xmin, Xmax])
    # ax.xaxis.set_major_locator(XMt)
    # ax.xaxis.set_minor_locator(Xmt)
    # ax.xaxis.set_major_formatter(X_format)
    ax.set_xlabel(r'Convergence time', rotation = 0, labelpad = 15, fontsize = 70)

    # Set frequency axes #
    # frq_min = 0.0
    # frq_max = 30000
    # frq_Mt = MultipleLocator(abs(frq_max - frq_min) / 5.00)  # frq major ticks
    # frq_mt = MultipleLocator(abs(frq_max - frq_min) / 10.0)  # frq minor ticks 
    # frq_format = FormatStrFormatter('%.2f')                  # integer labels
    # ax.set_ylim([frq_min, frq_max])
    # ax.yaxis.set_major_locator(frq_Mt)
    # ax.yaxis.set_minor_locator(frq_mt)
    # ax.yaxis.set_major_formatter(frq_format)
    ax.set_ylabel(r'Frequency', rotation = 90, labelpad = 15, fontsize = 70)

    # Set global features #
    plt.tick_params(which = 'major', direction = 'in', length = 18, width = 2.8, 
                    bottom = True, top = True, left = True, right = True, pad = 20)
    plt.tick_params(which = 'minor', direction = 'in', length = 12, width = 2.2, 
                    bottom = True, top = True, left = True, right = True)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(4.0)

    # Plot histogram #
    plt.hist(data, bins = Nbins, histtype = 'bar', align = 'left',
                   color = 'dodgerblue', density = False, zorder = 1)

    plt.title(par_value, x = 0.5, y = 1.01, fontsize = 70)
    plt.savefig(output_name, bbox_inches='tight')
    plt.clf()
    plt.close()

def __main():
    bins = 150
    input_name = get_input()
    data = data_parse(input_name)
    output_name = get_output()
    par_value = get_parameter()
    plot_histogram(data, bins, output_name, par_value)

if __name__ == "__main__":
    __main()