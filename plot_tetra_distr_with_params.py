import uproot
import matplotlib.pyplot as plt
import mplhep as hep
import pandas as pd
import numpy as np
import argparse
import contextlib
import subprocess
import Template_helper_methods

plt.style.use(hep.style.ROOT)
import matplotlib as mpl
mpl.rcParams['axes.labelsize'] = 40
mpl.rcParams['xaxis.labellocation'] = 'center'
mpl.rcParams['lines.linewidth'] = 2

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('commandFile', type=str)
    parser.add_argument('-s', '--suffix', default="")
    already_known = parser.add_argument_group("already known group", 
                                              "If you already know the parameters just use -p").add_mutually_exclusive_group(required=True)
    already_known.add_argument('-f', '--files', nargs='+', default=[])
    already_known.add_argument('-p', '--params', nargs=5, default=[], type=Template_helper_methods.CombineParam)
    args = parser.parse_args()
    
    zero_positions = np.zeros((len(args.files), 5)
                              ,dtype=float)
    
    positions = {
        "N":0,
        "RBW1":1,
        "RBW3":2,
        "RPhi12":3,
        "RPhi23":4
    }
    
    plot_command = [
        "python3",
        "create_1D_mass_interf_template_3_reso.py",
        "-c /eos/home-m/msrivast/www/data/Mass_Template_BasePlates/ghz1/CrossSections.csv",
        "-b ../Tetraquark_Analysis/Background/JPsi_Results/*_MELA_*.root",
        "-o /afs/cern.ch/user/m/msrivast/private/CMSSW_10_2_5/src/HiggsAnalysis/CombinedLimit",
        "-a 468 486 153",
        "-ba 824 5005"
    ]
    
    plot_command = " ".join(plot_command)
    print(plot_command)
    
    if args.files:
        with contextlib.ExitStack() as stack: #ExitStack takes in a bunch of things and closes them after so you don't have to!
            files = [
                stack.enter_context(uproot.open(fname)) for fname in args.files
            ]
            file_data = [file['limit'].arrays(list(positions.keys()) + ['deltaNLL'], library='pd') for file in files]
            
            for nf, file in enumerate(file_data):
                temp = file.loc[file['deltaNLL'].idxmin()][list(positions.keys())]
                # print(file.loc[file['deltaNLL'].idxmin()])
                runstr = ""
                for ne, element in enumerate(positions.keys()):
                    zero_positions[nf][positions[element]] = temp[element]
                    runstr += element + '=' + str(zero_positions[nf][positions[element]]) + " "
                    print(runstr)
                subprocess.run("python3 $(cat " + args.commandFile + ") -os " + runstr, shell=True)
                subprocess.run("mv scaled_hist.png local_files/scaled_hist_"+args.files[nf].split('/')[-1].split('.')[0]+args.suffix+'.png', shell=True)
    else:
        runstr = ""
        for param in args.params:
            runstr += param + " "
        print(runstr)
        subprocess.run("python3 $(cat " + args.commandFile + ") -os " + runstr, shell=True)
        params_in_name = sorted(args.params)
        subprocess.run("mv scaled_hist.png local_files/scaled_hist_" + "_".join(params_in_name) + args.suffix + '.png', shell=True)