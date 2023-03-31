import uproot
import matplotlib.pyplot as plt
import mplhep as hep
import pandas as pd
import numpy as np
import argparse
import contextlib
import Template_helper_methods as thm

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
    parser.add_argument('--log', action="store_true")
    parser.add_argument('-pre', '--prefix', default="")
    parser.add_argument('-nk', '--nokill', action="store_true")
    parser.add_argument('-max', '--maximum', type=float)
    args = parser.parse_args()
    
    zero_positions = np.zeros((len(args.files), len(args.branches))
                              ,dtype=float)
    
    with contextlib.ExitStack() as stack: #ExitStack takes in a bunch of things and closes them after so you don't have to!
        files = [
            stack.enter_context(uproot.open(fname)) for fname in args.files
        ]
        
        file_data = [file['limit'].arrays(args.branches + ['deltaNLL'], library='pd') for file in files]
        
        for ne, element in enumerate(args.branches):
            plt.cla()
            for nf, (data, name) in enumerate(zip(file_data, args.files)):
                data.sort_values(element, inplace=True, ignore_index=True)
                zero_positions[nf][ne] = data.loc[data['deltaNLL'].idxmin()][element]
                
                x = data[element]
                y=2*data['deltaNLL']
                if not args.nokill:
                    x, y = thm.killPoints(x,y)
                labelstr = "{:.3e}".format(zero_positions[nf][ne]) + '\n' + name.split('/')[-1].split('.')[0].split('_')[1].replace('tetra', '')
                plt.plot(x, y, lw=2, label=labelstr, linestyle='dashed')
                
            plt.gca().axhline(color='black', lw=2)
            plt.gca().axvline(color='black', lw=2)
            plt.grid(True)
            plt.ylabel(r'$2\times\Delta NLL$')
            plt.xlabel(element)
            plt.legend(loc='upper right')
            if args.maximum:
                plt.ylim(top=args.maximum)
            if args.log:
                plt.yscale('symlog')
                plt.ylim(bottom=0)
            plt.tight_layout()
            plt.savefig('local_files/' + args.prefix + element + '.png')
    
    # for file_params in zero_positions:
    #     file_params = list(map(str, file_params))
    #     passed_dict = {}
    #     for n, param in enumerate(args.branches):
    #         passed_dict[param] = file_params[n]
        
        # subprocess.run("python3 $(cat command.txt) -os " + file_params, shell=True)