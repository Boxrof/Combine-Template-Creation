import re
import tqdm
import glob
import uproot
import argparse
import numpy as np
import mplhep as hep
import Template_creator
import Template_helper_methods
import matplotlib.pyplot as plt
from collections.abc import Iterable

def place_that_list(filename):
    """Places different filenames according to their name when input. 
    Filenames should be either of form:
    BW<Reso Num>_pure if a pure sample
    BW<Reso Num>_phi_<phase>_BW<Reso Num>_phi_<phase>
    where phase is either 0 or pi_over_2

    Parameters
    ----------
    filename : string
        the filename that is being given
    """
    global insertionList
    if '_pure' in filename:
        if 'BW1' in filename:
            insertionList[0] = filename
        elif 'BW2' in filename:
            insertionList[1] = filename
        else: #BW3
            insertionList[2] = filename
    elif '_phi_pi_over_2' in filename:
        if 'BW2' in filename and 'BW1' in filename:
            insertionList[4] = filename
        elif 'BW2' in filename and 'BW3' in filename:
            insertionList[7] = filename
        else: #BW1 and BW3
            insertionList[6] = filename
    elif "_phi_0_" in filename and "_phi_0" in filename[filename.find('_phi_0')+6:]:
        if 'BW2' in filename and 'BW1' in filename:
            insertionList[3] = filename
        elif 'BW2' in filename and 'BW3' in filename:
            insertionList[8] = filename
        else: #BW1 and BW3
            insertionList[5] = filename
    else:
        print("whoops!")
        
def flatten(lst):
    for i in lst:
        if isinstance(i, Iterable) and not isinstance(i, (str, bytes)):
            yield from flatten(i)
        else:
            yield i

def UNIX_expansion_input(s):
    sPaths = glob.glob(s)
    if sPaths:
        return sPaths
    
    raise argparse.ArgumentError("File: " + s + " Not Found!")
        
def main(raw_args=None):
    global insertionList
    plt.style.use(hep.style.ROOT)

    parser = argparse.ArgumentParser()
    # parser.add_argument('filename')
    parser.add_argument('-n', '--nbins', default=120, type=int,
                        help="The number of bins you want")
    parser.add_argument('-o', '--outFolder', default='./',
                        help="The directory you'd like to output to")
    parser.add_argument('-na', '--name', default="Mass_Template", type=str,
                        help="What you'd like to name your folder")
    parser.add_argument('-c', '--crossSection', required=True,
                        help="The cross section file with all your sample names inside")
    parser.add_argument('-b', '--backgrounds', nargs='+', required=True, type=UNIX_expansion_input,
                        help="Your ROOT files containing the background samples")
    parser.add_argument('-ba', '--bkgAreas', nargs='+', required=True, type=float,
                        help="The areas for each background")
    parser.add_argument('-i', '--interference', action='store_true',
                        help="Turn this on if you are fitting for interference")
    scale_or_test = parser.add_mutually_exclusive_group()
    scale_or_test.add_argument('-os', '--outScaled', nargs='+', default=[], type=Template_helper_methods.CombineParam,
                        help="If this parameter is enabled ignore everything and plot the histogram based on the parameters given")
    scale_or_test.add_argument('-t', '--test', action="store_true",
                               help="Enable this option if you want to test your formulas")
    scale_or_test.add_argument('-a', '--animate', default=['', ''], nargs=2,
                               help="Enter a scan file from Combine and a variable to scan over to show an animation of the distribution as a function of the scan's progression ")
    args = parser.parse_args(raw_args)
    
    """
    Cross section files should be arranged as follows:
    <absolute file path>, <cross section>, <uncertainty in cross section> (last one is optional)
    The easiest to way to generate these is with the get_cross_section_from_LHE_file function in the lhe2root repo:
    https://github.com/hexutils/lhe2root/blob/main/lhe2root_methods.py
    """
    args.backgrounds = list(flatten(args.backgrounds))
    if len(args.bkgAreas) != len(args.backgrounds):
        print()
        print(args.bkgAreas)
        print(args.backgrounds)
        print()
        raise argparse.ArgumentError("Background argument and bkgArea arguments should be of the same length!")

    coupling_hunter = re.compile(r'\w+_ghzpzp(\d)_?\S+')

    data_samples = {}
    cross_section_samples = {}
    bkg_samples = {}
    with open(args.crossSection) as f:
        test = f.readline()
        for line in tqdm.tqdm(f):
            line = line.strip().split(',')
            with uproot.open(line[0]) as dataFile:
                dataFile = dataFile[dataFile.keys()[0]]
                sample = dataFile["M4L"].array(library="np")
                
                data_samples[line[0].split('/')[-1]] = sample
                cross_section_samples[line[0].split('/')[-1]] = float(line[1])
    
    insertionList = [None]*len(data_samples)
    
    for sampleName in data_samples.keys():
        place_that_list(sampleName)
    
    for bkg in args.backgrounds:
        with uproot.open(bkg) as dataFile:
            dataFile = dataFile[dataFile.keys()[0]]
            sample = dataFile["M4L"].array(library="np")
            bkg_samples[bkg.split('/')[-1].split('.')[0].split('_')[-1]] = sample
    
    
    Three_BW_Creation = None
    
    if args.interference:
        Three_BW_Creation = Template_creator.Interf_Reso_template_creator_1D(args.outFolder, args.name,
            bkg_samples.values(), bkg_samples.keys(), args.bkgAreas, 6, 9,
            *list(map(data_samples.get,insertionList)),
            *list(map(cross_section_samples.get, insertionList)),
            args.nbins)
    else:
        Three_BW_Creation = Template_creator.Reso_template_creator_1D(args.outFolder, args.name,
            bkg_samples.values(), bkg_samples.keys(), args.bkgAreas, 6, 9,
            *list(map(data_samples.get, insertionList[:3])),
            *list(map(cross_section_samples.get, insertionList[:3])), #only index to the third value so that you only get the pure terms
            args.nbins
        )
    
    if args.outScaled:
        params = {}
        for param in args.outScaled:
            param = param.split('=')
            params[param[0].lower()] = param[1]
        
        if len(args.outScaled) == 3 and not args.interference:
            Three_BW_Creation.histo_based_on_params(params['n1'], params['n2'], params['n3'])
        elif len(args.outScaled == 5) and args.interference:
            Three_BW_Creation.histo_based_on_params(params['n'], params['rbw1'], params['rbw3'], params['rphi12'], params['rphi23'])
        else:
            raise argparse.ArgumentError("--outScaled should have 5 arguments for the interference case, and 3 for the non-interference case!")
    elif args.test:
        print("testing...")
        Three_BW_Creation.check_for_correct_formulation()
    elif args.animate[0]:
        Three_BW_Creation.animate_over_scan(*args.animate)
    else:
        Three_BW_Creation.dump()
        Three_BW_Creation.create_datacards()
        Three_BW_Creation.stackPlot(args.nbins)
        # Three_BW_Creation.plot_overall_interference()
    
if __name__ == "__main__":
    main()
