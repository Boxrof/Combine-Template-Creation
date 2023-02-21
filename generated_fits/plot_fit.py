import uproot
import matplotlib.pyplot as plt
import mplhep as hep
import pandas as pd
import argparse
import contextlib

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
    parser.add_argument('-c', '--clean', action="store_false",
                        help="Activate this if you want all your plots on the same canvas")
    args = parser.parse_args()
    
    
    with contextlib.ExitStack() as stack: #ExitStack takes in a bunch of things and closes them after so you don't have to!
        files = [
            stack.enter_context(uproot.open(fname)) for fname in args.files
        ]
        
        file_data = [file['limit'].arrays(args.branches + ['deltaNLL'], library='pd') for file in files]
        
        for element in args.branches:
            if args.clean:
                plt.cla()
            for data, name in zip(file_data, args.files):
                data.sort_values(element, inplace=True, ignore_index=True)
                labelstr = name if args.clean else element
                plt.plot(data[element], 2*data['deltaNLL'], lw=2, label=labelstr)
                
            plt.gca().axhline(color='black', lw=2)
            plt.gca().axvline(color='black', lw=2)
            plt.grid(True)
            plt.ylabel(r'$2\times\Delta NLL$')
            if args.clean:
                plt.xlabel(element)
            plt.legend(loc='upper right')
            if args.clean:
                plt.tight_layout()
                plt.savefig(element + '.png')
                
        if not args.clean:
            plt.tight_layout()
            plt.savefig('Everything.png')
