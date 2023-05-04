import uproot
import matplotlib.pyplot as plt
import mplhep as hep
import pandas as pd
import numpy as np
import argparse
import contextlib
import Template_helper_methods as thm
from scipy import interpolate

plt.style.use(hep.style.ROOT)
import matplotlib as mpl
mpl.rcParams['axes.labelsize'] = 40
mpl.rcParams['xaxis.labellocation'] = 'center'
mpl.rcParams['lines.linewidth'] = 2

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('files', nargs='+')
    parser.add_argument('-b', '--branches', nargs='+', default=['RF'],
                        help="The branches you want to plot")
    parser.add_argument('--log', action="store_true",
                        help="If enabled, make the y axis log-linear")
    parser.add_argument('-pre', '--prefix', default="",
                        help="The prefix attached to the plots outputted")
    parser.add_argument('-nk', '--nokill', action="store_true",
                        help="If enabled, do not kill points")
    parser.add_argument('-max', '--maximum', type=float,
                        help="The maximum y value of the plots generated")
    parser.add_argument('-r', '--range', nargs=2, type=float,
                        help="The xrange of the plots generated in 2 arguments: <lower> <upper>")
    parser.add_argument('-i', '--interference', action="store_true",
                        help="Use this flag if you are looking at an interference sample")
    args = parser.parse_args()
    
    zero_positions = np.zeros((len(args.files), len(args.branches))
                              ,dtype=float)
    
    if args.range and (args.range[0] >= args.range[1]):
        raise argparse.ArgumentError("Ranges should go <lowest> <highest>")
    
    with contextlib.ExitStack() as stack: #ExitStack takes in a bunch of things and closes them after so you don't have to!
        files = [
            stack.enter_context(uproot.open(fname)) for fname in args.files
        ]
        
        file_data = [file['limit'].arrays(args.branches + ['deltaNLL'], library='pd') for file in files]
        
        
        for ne, element in enumerate(args.branches):
            plt.cla()
            for nf, (data, name) in enumerate(zip(file_data, args.files)):
                data.sort_values(element, inplace=True, ignore_index=True)
                
                central_value_index = data['deltaNLL'].idxmin()
                zero_positions[nf][ne] = data.loc[central_value_index][element]
                x = np.array(data[element])
                y = np.array(2*data['deltaNLL'])
                
                below_zero_point, correspondingx = y[:central_value_index], x[:central_value_index]
                below_zero_under = np.argmin(np.abs(np.where(below_zero_point - 1 < 0, np.inf, below_zero_point)))
                below_zero_over = np.argmin(np.abs(np.where(below_zero_point - 1 > 0, np.inf, below_zero_point)))
                lower_bound = thm.interpolate_uncertainty(1, 
                                                          correspondingx[below_zero_over], below_zero_point[below_zero_over],
                                                          correspondingx[below_zero_under], below_zero_point[below_zero_under])
                
                above_zero_point, correspondingx = y[central_value_index:], x[central_value_index:]
                above_zero_under = np.argmin(np.abs(np.where(above_zero_point - 1 < 0, np.inf, above_zero_point)))
                above_zero_over = np.argmin(np.abs(np.where(above_zero_point - 1 > 0, np.inf, above_zero_point)))
                upper_bound = thm.interpolate_uncertainty(1, 
                                                          correspondingx[below_zero_under], above_zero_point[below_zero_under],
                                                          correspondingx[below_zero_over], above_zero_point[below_zero_over])
                
                index_above = np.argmin(np.abs(above_zero_point - 1)) #finds the point that is closest to 1 for the error bound
                
                plt.gca().axvline(lower_bound, linestyle='dashed', color='black')
                plt.gca().axvline(upper_bound, linestyle='dashed', color='black')
                
                lower_uncertainty = np.abs(zero_positions[nf][ne] - lower_bound)
                upper_uncertainty = np.abs(zero_positions[nf][ne] - upper_bound) 
                
                if not args.nokill:
                    x, y = thm.killPoints(x,y)
                
                labelstr = r"${{{:.2f}}}^{{+{:.2f}}}_{{-{:.2f}}}$".format(
                    zero_positions[nf][ne], upper_uncertainty, lower_uncertainty)
                plt.plot(x, y, lw=2, label=labelstr, linestyle='dashed')
                
            plt.gca().axhline(color='black', lw=2)
            plt.gca().axvline(color='black', lw=2)
            plt.grid(True)
            plt.ylabel(r'$2\times\Delta NLL$')
            plt.xlabel(element)
            plt.legend(loc='best')
            if args.maximum:
                plt.ylim(top=args.maximum)
            if args.log:
                plt.yscale('symlog')
                plt.ylim(bottom=0)
            if args.range:
                plt.xlim(args.range[0], args.range[1])
            plt.tight_layout()
            
            placement = ""
            if args.interference:
                placement = 'local_files/interference/'
            else:
                placement = 'local_files/non_interference/'
            plt.savefig(placement + args.prefix + element + '.png')
    
    # for file_params in zero_positions:
    #     file_params = list(map(str, file_params))
    #     passed_dict = {}
    #     for n, param in enumerate(args.branches):
    #         passed_dict[param] = file_params[n]
        
        # subprocess.run("python3 $(cat command.txt) -os " + file_params, shell=True)