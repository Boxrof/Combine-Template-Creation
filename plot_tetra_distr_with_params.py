import uproot
import matplotlib.pyplot as plt
import mplhep as hep
import pandas as pd #required for uproot output
import numpy as np
import argparse
import contextlib
import subprocess
import Template_helper_methods
import matplotlib as mpl
import create_1D_mass_interf_template_3_reso as reso_temp
import tqdm

def main(raw_args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('commandFile', type=str,
                        help="This is the file that contains the command for plotting")
    parser.add_argument('-s', '--suffix', default="")
    parser.add_argument('-n', '--nbins', default="")
    parser.add_argument('-i', '--interference', action='store_true')
    already_known = parser.add_argument_group("already known group", 
                                              "If you already know the parameters just use -p").add_mutually_exclusive_group(required=True)
    already_known.add_argument('-f', '--files', nargs='+', default=[])
    already_known.add_argument('-p', '--params', nargs=5, default=[], type=Template_helper_methods.CombineParam)
    args = parser.parse_args(raw_args)
    
    zero_positions = np.zeros((len(args.files), 5)
                              ,dtype=float)
    positions = {}
    
    if args.interference:
        positions = {
            "N":0,
            "RBW1":1,
            "RBW3":2,
            "RPhi12":3,
            "RPhi23":4
        }
    else:
        positions = {
            "N1":0,
            "N1":1,
            "N2":2,
            "N3":3
        }
    if args.files:
        with contextlib.ExitStack() as stack: #ExitStack takes in a bunch of things and closes them after so you don't have to!
            files = [
                stack.enter_context(uproot.open(fname)) for fname in args.files
            ]
            file_data = [file['limit'].arrays(list(positions.keys()) + ['deltaNLL'], library='pd') for file in files]
            
            for nf, file in enumerate(file_data):
                temp = file.loc[file['deltaNLL'].idxmin()][list(positions.keys())]
                param_vals = []
                for element in tqdm.tqdm(positions.keys()):
                    zero_positions[nf][positions[element]] = temp[element]
                    param_vals.append(element + '=' + str(zero_positions[nf][positions[element]]))
                # subprocess.run("python3 $(cat " + args.commandFile + ") -os " + runstr, shell=True)
                with open(args.commandFile) as f:
                    input_args = f.read().strip().split()
                    if args.nbins:
                        reso_temp.main(input_args + ["-os"] + param_vals + ["-n"] + [args.nbins])
                    else:
                        reso_temp.main(input_args + ["-os"] + param_vals)
                    subprocess.run("mv scaled_hist.png local_files/scaled_hist_"+args.files[nf].split('/')[-1].split('.')[0]+"_"+args.suffix+'.png', shell=True)
                    subprocess.run("mv scaled_hist_with_bkg.png local_files/scaled_hist_" + args.files[nf].split('/')[-1].split('.')[0]+"_" + "_" + args.suffix + '_with_bkg.png', shell=True)
    else:
        param_vals = []
        for param in args.params:
            param_vals.append(param)
        # subprocess.run("python3 $(cat " + args.commandFile + ") -os " + runstr, shell=True)
        with open(args.commandFile) as f:
            input_args = f.read().strip().split()
            if args.nbins:
                reso_temp.main(input_args + ["-os"] + param_vals + ["-n"] + [args.nbins])
            else:
                reso_temp.main(input_args + ["-os"] + param_vals)
            params_in_name = sorted(args.params)
            subprocess.run("mv scaled_hist.png local_files/scaled_hist_" + "_".join(params_in_name) + "_" + args.suffix + '.png', shell=True)
            subprocess.run("mv scaled_hist_with_bkg.png local_files/scaled_hist_" + "_".join(params_in_name) + "_" + args.suffix + '_with_bkg.png', shell=True)

if __name__ == "__main__":
    plt.style.use(hep.style.ROOT) #use this style initially
    mpl.rcParams['axes.labelsize'] = 40 #these overwrite the style set initially
    mpl.rcParams['xaxis.labellocation'] = 'center'
    mpl.rcParams['lines.linewidth'] = 2
    
    main()