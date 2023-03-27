import Mass_interference_helper_methods as mihm
from collections.abc import Iterable
import matplotlib.pyplot as plt
import Template_helper_methods
import mplhep as hep
import pandas as pd
import numpy as np
import contextlib
import warnings
import shutil
import uproot
import tqdm
import ROOT
import time
import copy
import os

plt.style.use(hep.style.ROOT)

class Template_creator(object):
    def __init__(self, output_directory, fname, bkgs, bkgNames, bkg_areas, lowerlim, upperlim):
        """This serves as a parent class for all other templates made

        Parameters
        ----------
        output_directory : str
            The directory you would like to output all your data in
        fname : str
            The filenames contained in all your outputs
        bkgs : list[list[float]]
            A list of background iterables. These can be lists of mass distriutions, lists of angle distributions, etc. bkgs is essentially a list of lists.
        bkgNames : list[str]
            A list of names for each background. This should be the same dimension as bkgs.
        bkg_areas : list[Union[float, int]]
            A list of areas for each background that you would like to scale to
        lowerlim : float
            The lower limit of your attribute's range (i.e. if it was phi, lowerlim would be -pi)
        upperlim : float
            The upper limit of your attribute's range (i.e. if it was phi, upperlim would be pi)
        """
        self.running_location = os.getcwd()
        if 'HexUtils' in self.running_location:
            self.HexUtils_path = self.running_location[:self.running_location.find('HexUtils') + 8] #8 is the number of characters in "HexUtils"
        else:
            self.HexUtils_path = ""
            warnings.warn("Not in an instance of HexUtils! This will lead to undefined behavior!")
        
        self.dimension = 0 #This attribute is set when creating a certain dimension of template
        
        
        self.output_directory = os.path.abspath(output_directory)
        if self.output_directory[-1] != '/':
            self.output_directory += '/' #make sure there is a slash at the end to ensure the shell knows it is a directory
        self.fname = fname.split('.')[0] #make sure there are no extensions in the filename
        
        self.lowerlim = lowerlim
        self.upperlim = upperlim
        self.discr_range = (0,1) #This is (0,1) for most templates, but (-1,1) for hypothesis interference templates
        
        self.signals = {} #this is a dictionary with values that look like (<raw_data>, <area desired>)
        self.scaled_signals = {} #this is a dictionary that should contain scaled histograms of your raw data.
        #so the values should look like (<counts>, <bins>)
        self.signal_weights = {} #this is a dictionary that should hold the probabilities from MELA should you choose to use them
        self.discr_signals = {} #this is a dictionary that should hold the discriminants that your signals are using
        #both signal weights and discr_weights contain ITERABLES! Each value should be an iterable of values!!
        
        ################### The dictionaries below are the versions of the above dictionaries but for background ###########
        self.bkgs = {}
        self.scaled_bkgs = {}
        self.bkg_weights = {}
        self.discr_bkgs = {}
        
        for name, bkg_sample, bkg_area in zip(bkgNames, bkgs, bkg_areas):
            self.bkgs[name] = (bkg_sample, bkg_area) #preprocessing the background samples
            
    def create_datacards(self, verbose=False, clean=True):
        """This function uses DatacardMaker_OnShell and MakeInputRoot_OnShell to generate the datacards and input for Higgs Combine

        Parameters
        ----------
        verbose : bool, optional
            Whether you would like verbosity, by default False
        clean : bool, optional
            Whether you would like to wipe your output folders before putting anything else in them, by default True
        
        Raises
        ------
        FileNotFoundError
            If there is no file, raises an error
        """
        filename = self.output_directory + self.fname + ".root"
        if not os.path.isfile(filename):
            raise FileNotFoundError("Make sure to dump the Template to a ROOT file first!")
        
        in_folder = self.output_directory + self.fname
        out_folder = in_folder+"_out"
        
        in_folder += '/'
        out_folder += '/'
        
        if not os.path.isdir(in_folder):
            os.mkdir(in_folder)
        if clean:
            os.system("rm " + in_folder + '*')
        shutil.move(filename, in_folder)
        
        if not os.path.isdir(out_folder):
            os.mkdir(out_folder)
        if clean:
            os.system("rm " + out_folder + '*')  
        
        runstr = "python3 MakeInputRoot_OnShell.py "
        runstr += in_folder + " " + out_folder
        if not verbose:
            runstr += " > /dev/null"
        os.system(runstr)
        # print(runstr)
        
        runstr = "python3 DatacardMaker_OnShell.py "
        runstr += out_folder
        if not verbose:
            runstr += "> /dev/null"
        os.system(runstr)
        # print(runstr)
    
    def scale_and_add_bkgs(self):
        """This is a placeholder for the same function in the 1D and 2D template cases
        """
        return 
    
    def stackPlot(self, nbins=40):
        """Generates a stack plot of your samples

        Parameters
        ----------
        nbins : int, optional
            The number of plots you want, by default 40
        """
        bins = []
        weights = []
        labels = []
        for bkg_name, (bkg, area) in self.bkgs.items():
            if any(bins):
                bkg, _ = np.histogram(bkg, bins=bins, range=(self.lowerlim, self.upperlim))
            else:
                bkg, bins = np.histogram(bkg, bins=nbins, range=(self.lowerlim, self.upperlim))
            bkg = Template_helper_methods.scale(area, bkg)
            labels.append(bkg_name)
            weights.append(bkg)
        
        for sig_name, (sig, area) in self.signals.items():
            sig, _ = np.histogram(sig, bins=bins, range=(self.lowerlim, self.upperlim))
            sig = Template_helper_methods.scale(area, bkg)
            plt.cla()
            hep.histplot(weights + [sig], bins=bins, label=labels + [sig_name], stack=True, lw=3)
            plt.legend()
            plt.savefig(self.output_directory + sig_name + '_stack.png')
        
class Template_Creator_1D(Template_creator):
    def __init__(self, output_directory, fname, bkgs, bkgNames, bkg_areas, lowerlim, upperlim):
        """This initialization takes in all the same inputs as the parent class.

        Parameters
        ----------
        output_directory : str
            The directory you would like to output all your data in
        fname : str
            The filenames contained in all your outputs
        bkgs : list[list[float]]
            A list of background iterables. These can be lists of mass distriutions, lists of angle distributions, etc. bkgs is essentially a list of lists.
        bkgNames : list[str]
            A list of names for each background. This should be the same dimension as bkgs.
        bkg_areas : list[Union[float, int]]
            A list of areas for each background that you would like to scale to
        lowerlim : float
            The lower limit of your attribute's range (i.e. if it was phi, lowerlim would be -pi)
        upperlim : float
            The upper limit of your attribute's range (i.e. if it was phi, upperlim would be pi)
        """
        super().__init__(output_directory, fname, bkgs, bkgNames, bkg_areas, lowerlim, upperlim)
        self.dimension = 1 #resets the dimension to 1
        
    def scale_and_add_bkgs(self, bins=40, scaleTo=True):
        """This is the 1-dimensional version of the function. 
        It serves to bin and scale the backgrounds given to their respective areas, then add them into one histogram.

        Parameters
        ----------
        bins : int, optional
            Either the number of bins you want, or a list of the bins you want, by default 40
        scaleTo : bool, optional
            If true, this function will scale the backgrounds before adding them, by default True

        Returns
        -------
        Tuple[numpy.ndarray, numpy.ndarray]
            an overall histogram pair of (counts, bins) a la a numpy histogram
        """
        names_samples_and_areas = list(self.bkgs.items())
            
        name, (sample, area) = names_samples_and_areas[0]
        bkg_sample, bins = np.histogram(sample, range=(self.lowerlim, self.upperlim), bins=bins)
        if scaleTo:
            bkg_sample = Template_helper_methods.scale(area, bkg_sample)
            self.scaled_bkgs[name] = (bkg_sample, bins)
        overall = bkg_sample
        
        for name, (sample, area) in names_samples_and_areas[1:]:
            bkg_sample, _ = np.histogram(sample, range=(self.lowerlim, self.upperlim), bins=bins)
            if scaleTo:
                bkg_sample = Template_helper_methods.scale(area, bkg_sample)
                self.scaled_bkgs[name] = (bkg_sample, bins)
            overall += bkg_sample
        
        return overall, bins

class Template_Creator_2D(Template_creator):
    def __init__(self, output_directory, fname, bkgs, bkgNames, bkg_areas, lowerlim, upperlim):
        """This initialization takes in all the same inputs as the parent class.

        Parameters
        ----------
        output_directory : str
            The directory you would like to output all your data in
        fname : str
            The filenames contained in all your outputs
        bkgs : list[list[float]]
            A list of background iterables. These can be lists of mass distriutions, lists of angle distributions, etc. bkgs is essentially a list of lists.
        bkgNames : list[str]
            A list of names for each background. This should be the same dimension as bkgs.
        bkg_areas : list[Union[float, int]]
            A list of areas for each background that you would like to scale to
        lowerlim : float
            The lower limit of your attribute's range (i.e. if it was phi, lowerlim would be -pi)
        upperlim : float
            The upper limit of your attribute's range (i.e. if it was phi, upperlim would be pi)
        """
        super().__init__(output_directory, fname, bkgs, bkgNames, bkg_areas, lowerlim, upperlim)
        self.dimension = 2   
        
    def scale_and_add_bkgs(self, bins=40, scaleTo=True):
        """This is the 2-dimensional version of the function. 
        It serves to bin and scale the backgrounds given to their respective areas, then add them into one histogram.

        Parameters
        ----------
        bins : Union[int, array_like], optional
            Either the number of bins you want, or a list of the bins you want, by default 40
        scaleTo : bool, optional
            If true, this function will scale the backgrounds before adding them, by default True

        Returns
        -------
        Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]
            an overall histogram pair of (counts, binsx, binsy) a la a numpy histogram
        """
        overall = np.ndarray((2,2))
        keys_to_follow = list(self.bkgs.keys())
        
        current_key = keys_to_follow[0]
        sample, area = self.bkgs[current_key]
        discr = self.discr_bkgs[current_key]
        
        bkg_sample, binsx, binsy = np.histogram2d(sample, discr, bins, range=[(self.lowerlim, self.upperlim), self.discr_range])
        if scaleTo:
            bkg_sample = Template_creator.scale(area, bkg_sample)
            self.scaled_bkgs[current_key] = (bkg_sample, binsx, binsy)
        overall = bkg_sample
        
        for current_key in keys_to_follow[1:]:
            sample, area = self.bkgs[current_key]
            discr = self.discr_bkgs[current_key]
            
            bkg_sample, _, _ = np.histogram2d(sample, discr, (binsx, binsy), range=[(self.lowerlim, self.upperlim), self.discr_range])
            if scaleTo:
                bkg_sample = Template_creator.scale(area, bkg_sample)
                self.scaled_bkgs[current_key] = (bkg_sample, binsx, binsy)
            overall += bkg_sample
            
        return overall, binsx, binsy 

class Interf_Coupling_template_creator(Template_Creator_2D): #WIP
    def __init__(self, output_directory, fname, bkgs, bkgNames, bkg_areas, lowerlim, upperlim,
                 pure1_weights, pure2_weights, interf_weights, weight_of_generation_hypothesis, mass_iterable,
                 interf_name, pure1_name, pure2_name,
                 bkg_pure1_weights, bkg_pure2_weights, bkg_interf_weights):
        """CURRENTLY A WORK IN PROGRESS. Designed to be a 2d template between different hypotheses

        Parameters
        ----------
        output_directory : str
            The directory you would like to output all your data in
        fname : str
            The filenames contained in all your outputs
        bkgs : list[list[float]]
            A list of background iterables. These can be lists of mass distriutions, lists of angle distributions, etc. bkgs is essentially a list of lists.
        bkgNames : list[str]
            A list of names for each background. This should be the same dimension as bkgs.
        bkg_areas : list[Union[float, int]]
            A list of areas for each background that you would like to scale to
        lowerlim : float
            The lower limit of your attribute's range (i.e. if it was phi, lowerlim would be -pi)
        upperlim : float
            The upper limit of your attribute's range (i.e. if it was phi, upperlim would be pi)
        pure1_weights : list[float]
            The weights for the first pure sample. Should be an iterable!
        pure2_weights : list[float]
            The weights for the second pure sample. Should be an iterable!
        interf_weights : list[float]
            The weights for the sample of the interference between the 2. Should be an iterable!
        weight_of_generation_hypothesis : str
            The weight for what the sample was generated at. This is used to normalize weights.
        mass_iterable : list[float]
            This is the mass iterable that you are using for one axis of the template
        interf_name : str
            The name of your interference sample
        pure1_name : str
            The name of your first pure sample
        pure2_name : str
            The name of your second pure sample
        bkg_pure1_weights : list[float]
            The weights of your background to your first pure sample's hypothesis
        bkg_pure2_weights : list[float]
            The weights of your background to your second pure sample's hypothesis
        bkg_interf_weights : list[float]
            The weights of your background to your interference sample's hypothesis
        """
        super().__init__(output_directory, fname, bkgs, bkgNames, bkg_areas, lowerlim, upperlim)
        
        self.discr_range = (-1,1)
        
        weight_of_generation_hypothesis = np.array(weight_of_generation_hypothesis)
        self.signal_weights[pure1_name] = np.array(pure1_weights)/weight_of_generation_hypothesis
        self.signal_weights[pure2_name] = np.array(pure2_weights)/weight_of_generation_hypothesis
        
        self.signal_weights[interf_name] = np.array(interf_weights)/weight_of_generation_hypothesis
        self.signal_weights[interf_name] -= self.weights[pure1_name] + self.weights[pure2_name]
        
        mass_iterable = np.array(mass_iterable)
        self.signals.update(dict.fromkeys([interf_name, pure1_name, pure2_name], mass_iterable ))
        
        d_interference_signal = self.signal_weights[interf_name]
        d_interference_signal /= 2*np.sqrt(self.signal_weights[pure1_name]*self.signal_weights[pure2_name])
        
        d_interference_bkg = bkg_interf_weights
        d_interference_bkg /= 2*np.sqrt(bkg_pure1_weights*bkg_pure2_weights)
        
class Interf_Reso_template_creator_1D(Template_Creator_1D):
    def __init__(self, output_directory, fname, bkgs, bkgNames, bkg_areas, lowerlim, upperlim,
                 BW1_0_0, BW2_0_0, BW3_0_0, BW12_0_0, BW12_05_0, BW13_0_0, BW13_0_05, BW23_0_0, BW23_0_05,
                 CS_BW1, CS_BW2, CS_BW3, CS_BW12_0_0, CS_BW12_05_0, CS_BW13_0_0, CS_BW13_0_05, CS_BW23_0_0, CS_BW23_0_05,
                 nbins, area1, area2, area3):
        """Initializes the 1D template for mass interference

        Parameters
        ----------
        output_directory : str
            The directory you would like to output all your data in
        fname : str
            The filenames contained in all your outputs
        bkgs : list[list[float]]
            A list of background iterables. These can be lists of mass distriutions, lists of angle distributions, etc. bkgs is essentially a list of lists.
        bkgNames : list[str]
            A list of names for each background. This should be the same dimension as bkgs.
        bkg_areas : list[Union[float, int]]
            A list of areas for each background that you would like to scale to
        lowerlim : float
            The lower limit of your attribute's range (i.e. if it was phi, lowerlim would be -pi)
        upperlim : float
            The upper limit of your attribute's range (i.e. if it was phi, upperlim would be pi)
        BW1_0_0 : list[float]
            The first pure sample. An iterable of masses.
        BW2_0_0 : list[float]
            The second pure sample. An iterable of masses.
        BW3_0_0 : list[float]
            The third pure sample. An iterable of masses.
        BW12_0_0 : list[float]
            The interference sample with 0 phase between 1 & 2. An iterable of masses.
        BW12_05_0 : list[float]
            The interference sample with pi/2 phase between 1 & 2 (hence the "05" indicating 0.5pi). An iterable of masses.
        BW13_0_0 : list[float]
            The interference sample with 0 phase between 1 & 3. An iterable of masses.
        BW13_0_05 : list[float]
            The interference sample with pi/2 phase between 1 & 3. An iterable of masses.
        BW23_0_0 : list[float]
            The interference sample with 0 phase between 2 & 3. An iterable of masses.
        BW23_0_05 : list[float]
            The interference sample with pi/2 phase between 2 & 3. An iterable of masses.
        CS_BW1 : float
            The cross section of the first pure sample
        CS_BW2 : float
            The cross section of the second pure sample
        CS_BW3 : float
            The cross section of the third pure sample
        CS_BW12_0_0 : float
            The cross section of the variable's name
        CS_BW12_05_0 : float
            The cross section of the variable's name
        CS_BW13_0_0 : float
            The cross section of the variable's name
        CS_BW13_0_05 : float
            The cross section of the variable's name
        CS_BW23_0_0 : float
            The cross section of the variable's name
        CS_BW23_0_05 : float
            The cross section of the variable's name
        nbins : int
            The cross section of the variable's name
        area1 : float
            The cross section of the variable's name
        area2 : float
            The cross section of the variable's name
        area3 : float
            The cross section of the variable's name
        """
        super().__init__(output_directory, fname, bkgs, bkgNames, bkg_areas, lowerlim, upperlim)
        self.string_forms = ["BW1", "BW2", "BW3", 
                        "BW1BW2_0_0", "BW1BW2_0.5_0", "BW1BW3_0_0", "BW1BW3_0_0.5", "BW2BW3_0_0", "BW2BW3_0_0.5"]
        #The 0.5 is for the physics model naming scheme
        
        self.final_scaling_funcs = {#this is how the scaling goes for all 5 parameters
            "BW1" : (lambda N, f1, f3, phi12, phi23: N*f1), 
            "BW2" : (lambda N, f1, f3, phi12, phi23: N*(1 - f1 - f3)), 
            "BW3" : (lambda N, f1, f3, phi12, phi23: N*f3), 
            "BW1BW2_0_0" : (lambda N, f1, f3, phi12, phi23 : N*np.sqrt(f1*(1-f1-f3))*np.cos(phi12)), 
            "BW1BW2_0.5_0" : (lambda N, f1, f3, phi12, phi23 : N*np.sqrt(f1*(1-f1-f3))*np.sin(phi12)), 
            "BW2BW3_0_0" : (lambda N, f1, f3, phi12, phi23 : N*np.sqrt((1-f1-f3)*f3)*np.cos(phi23)), 
            "BW2BW3_0_0.5" : (lambda N, f1, f3, phi12, phi23 : N*np.sqrt((1-f1-f3)*f3)*np.sin(phi23)), 
            "BW1BW3_0_0" : (lambda N, f1, f3, phi12, phi23 : N*np.sqrt(f1*f3)*np.cos((phi12 - phi23))), 
            "BW1BW3_0_0.5" : (lambda N, f1, f3, phi12, phi23 : N*np.sqrt(f1*f3)*np.sin((phi23 - phi12)))
            }
        
        normalization_factors = np.zeros(3, dtype=float)
        
        
        self.signals["BW1"] = (np.array(BW1_0_0), CS_BW1)
            
        BW1_0_0, bins = np.histogram(BW1_0_0, bins=nbins, range=(lowerlim, upperlim))
        BW1_0_0 = Template_helper_methods.scale(CS_BW1, BW1_0_0)
        # BW1_0_0, normalization_factors[0] = Template_helper_methods.scale(1, BW1_0_0, return_scale_factor=True)
        # BW1_0_0 /= CS_BW1
        
        self.scaled_signals["BW1"] = (BW1_0_0/CS_BW1, bins)
        
        self.signals["BW2"] = (np.array(BW2_0_0), CS_BW2)
        BW2_0_0, _ = np.histogram(BW2_0_0, bins=bins, range=(lowerlim, upperlim))
        BW2_0_0 = Template_helper_methods.scale(CS_BW2, BW2_0_0)
        # BW2_0_0, normalization_factors[1] = Template_helper_methods.scale(1, BW2_0_0, return_scale_factor=True)
        # BW2_0_0 /= CS_BW2
        
        self.scaled_signals["BW2"] = (BW2_0_0/CS_BW2, bins)
        
        self.signals["BW3"] = (np.array(BW3_0_0), CS_BW3)
        BW3_0_0, _ = np.histogram(BW3_0_0, bins=bins, range=(lowerlim, upperlim))
        BW3_0_0 = Template_helper_methods.scale(CS_BW3, BW3_0_0)
        # BW3_0_0, normalization_factors[2] = Template_helper_methods.scale(1, BW3_0_0, return_scale_factor=True)
        # BW3_0_0 /= CS_BW3
        
        self.scaled_signals["BW3"] = (BW3_0_0/CS_BW3, bins)
        
        interfList = [BW12_0_0, BW12_05_0, BW13_0_0, BW13_0_05, BW23_0_0, BW23_0_05] #list of all the interference terms
        interfCSList = [CS_BW12_0_0, CS_BW12_05_0, CS_BW13_0_0, CS_BW13_0_05, CS_BW23_0_0, CS_BW23_0_05]
        
        print(np.sum(BW1_0_0), np.sum(BW2_0_0), np.sum(BW3_0_0), normalization_factors)
        for n, interference_term in enumerate(interfList):
            
            interference_term, _ = np.histogram(interference_term, bins=bins, range=(lowerlim, upperlim))            
            self.signals[self.string_forms[n+3]] = (np.array(interference_term), interfCSList[n])
            interference_term = Template_helper_methods.scale(interfCSList[n], interference_term)
            
            dividing_to_normalize = 1
            subtraction_terms = np.zeros(len(BW1_0_0))
            
            if "BW1BW2" in self.string_forms[n + 3]:
                dividing_to_normalize = np.sqrt(CS_BW1*CS_BW2)
                print("12:", dividing_to_normalize)
                subtraction_terms += BW1_0_0 + BW2_0_0
            elif "BW1BW3" in self.string_forms[n + 3]:
                dividing_to_normalize = np.sqrt(CS_BW1*CS_BW3)
                print("13:", dividing_to_normalize)
                subtraction_terms += BW1_0_0 + BW3_0_0
            elif "BW2BW3" in self.string_forms[n+3]:
                dividing_to_normalize = np.sqrt(CS_BW2*CS_BW3)
                print("23:", dividing_to_normalize)
                subtraction_terms += BW2_0_0 + BW3_0_0
            else:
                print("ERROR")
                exit()
            
            print("interference area of", "{:.2e}".format(np.sum(interference_term)), "for", self.string_forms[n + 3])
            interference_term -= subtraction_terms
            print("Interference area subtracted down to", "{:.2e}".format(np.sum(interference_term)))
            interference_term = interference_term/dividing_to_normalize
            print("Final area of interference: ", np.sum(interference_term), np.sum(np.abs(interference_term)), '\n')
            self.scaled_signals[self.string_forms[n+3]] = (interference_term, bins)
            
        self.bins = bins
        
    def dump(self):
        """Dumps the created items in a file at self.output_directory/self.fname
        """
        with uproot.recreate(self.output_directory + self.fname + ".root") as f:
            
            for signal in self.scaled_signals.keys():
                
                if signal == "BW1" or signal == "BW2" or signal == "BW3":
                    f["ggH_0PM_"+signal] = self.scaled_signals[signal]
                        
                # # BW2_0_0 = Template_helper_methods.scale(BW2_0_0, area2)
                # f["ggH_0PM_BW2"] = self.scaled_signals["BW2"]
            
                # # BW3_0_0 = Template_helper_methods.scale(BW3_0_0, area3)
                # f["ggH_0PM_BW3"] = self.scaled_signals["BW3"]
                
                else:
                    interference_term, bins = self.scaled_signals[signal]
                    
                    pos = np.maximum(interference_term.copy(),0) #splits the template up into positive and negative as you're supposed to
                    neg = -1*np.minimum(interference_term.copy(),0)
                    
                    if np.any(pos):
                        f["ggH_0PM_" + signal + "_positive"] = (pos, bins)
                    if np.any(neg):
                        f["ggH_0PM_" + signal + "_negative"] = (neg, bins)

            f["bkg_ggzz"] = self.scale_and_add_bkgs(bins, scaleTo=True)
    
    def histo_based_on_params(self, N, f1, f3, phi12, phi23, fname="scaled_hist"):
        """Returns a histogram based on the 5 parameters that are fitted in the template
        While not necessarily part of creating the template, this is a useful place to put this function

        Parameters
        ----------
        f1 : float
            The fraction for BW1
        f2 : float
            The fraction for BW2
        f3 : float
            The fraction for BW3
        phi12 : float
            The phase between BW1 and BW2
        phi23 : float
            The phase between BW2 and BW3
        fname : str
            The filename you want, by default "scaled_hist"

        Returns
        -------
        Tuple[list[float], list[float]]
            A numpy histogram of everything
        """
        param_dict = locals() #this stores the parameter names
        
        params = [N, f1, f3, phi12, phi23]
        total = np.zeros(len(self.bins) - 1, dtype=float)
        params = list(map(float, params))
        plt.figure()
        for signal in self.scaled_signals.keys():
            temp_counts = self.final_scaling_funcs[signal](*params)*self.scaled_signals[signal][0]
            # temp_counts, _ = Template_helper_methods.scale(scalefactor, *self.scaled_signals[signal])
            hep.histplot(temp_counts, self.bins, label=signal, lw=3)
            total += temp_counts
            
        if np.any(total < 0):
            warnings.warn("Bins Cannot be Negative! Fix your Methodology!", RuntimeWarning)
        
        hep.histplot(total, self.bins, lw=4, label="all:" + "{:.1f}".format(np.sum(total)), color="black")
        plt.gca().axhline(lw=2, color='black')
        plt.legend()
        titlestr = ""
        for name, i in param_dict.items():
            try:
                if "phi" in name:
                    titlestr += name + ": {:.2f} ".format(float(i))# + r"$\pi$ "
                else:
                    titlestr += name + ": {:.2f} ".format(float(i))
            except:
                pass
        plt.title(titlestr)
        plt.xlabel(r"$m_{4\mu}[GeV]$")
        plt.savefig(fname+".png")
        return total, self.bins
    
    def check_for_correct_formulation(self, iters=10):
        """Checks Whether your formulation EVER returns a negative bin!

        Parameters
        ----------
        iters : int, optional
            The number of times you sample the range of events, by default 10
        Returns
        -------
        bool
            Whether your functions work!
        """
        os.system("rm failure*.png")
        
        allWorked = True
        failCount = 1
        
        for f1 in tqdm.tqdm(np.linspace(1, 0, iters, endpoint=False), desc="f1"):
            for f3 in tqdm.tqdm(np.linspace(1 - f1, 0, iters, endpoint=False), desc="f3", leave=False):
                for phi12 in tqdm.tqdm(np.linspace(-np.pi, np.pi, iters), desc="phi1", leave=False):
                    for phi23 in tqdm.tqdm(np.linspace(-np.pi, np.pi, iters), desc="phi1", leave=False):
                        total = np.zeros(len(self.bins) - 1, dtype=float)
                        params = [2000, f1, f3, phi12, phi23]
                        for signal in self.scaled_signals.keys():
                            total += self.final_scaling_funcs[signal](*params)*self.scaled_signals[signal][0]
                        if np.any(total < -1):
                            print("failed at f1=",f1,"f3=",f3,"phi1=",phi12,"phi3=",phi23)
                            self.histo_based_on_params(*params, fname="failure" + str(failCount))
                            failCount += 1
                            allWorked = False
        return allWorked
        
    def plot_overall_interference(self):
        """This plots all the different combinations of the three phases
        """
        # print(self.scaled_signals)
        pures = ([key for key in self.scaled_signals.keys() if ("1BW" not in key and "1BW" not in key and "2BW" not in key)], 
                 [value for key, value in self.scaled_signals.items() if ("1BW" not in key and "1BW" not in key and "2BW" not in key)])
        for interf12 in tqdm.tqdm(["BW1BW2_0_0", "BW1BW2_0.5_0"], desc="Top Level of interference loop"):
            for interf13 in tqdm.tqdm(["BW1BW3_0_0", "BW1BW3_0_0.5"], leave=False, desc="Second Level of interference loop"):
                for interf23 in tqdm.tqdm(["BW2BW3_0_0", "BW2BW3_0_0.5"], leave=False, desc="Bottom Level of interference loop"):
                    names, terms = copy.deepcopy(pures)
                    names += [interf12, interf13, interf23]
                    # print(names)
                    terms += [self.scaled_signals[interf12], self.scaled_signals[interf13], self.scaled_signals[interf23]]

                    mihm.plot_overall_interference(terms, names, 
                                                   self.output_directory, interf12+"_"+interf13+"_"+interf23+"_"+self.fname)

class Significance_Hypothesis_template_creator_1D(Template_Creator_1D):
    def __init__(self, output_directory, fname, bkgs, bkgNames, bkg_areas, lowerlim, upperlim,
                 signal1, signal1_name, signal2, signal2_name, signal_area, nbins):
        super().__init__(output_directory, fname, bkgs, bkgNames, bkg_areas, lowerlim, upperlim)

        with uproot.recreate(self.output_directory + self.fname + ".root") as f:
            self.signals[signal1_name] = signal1
            signal1, bins = np.histogram(signal1, bins=nbins, range=(lowerlim, upperlim))
            signal1 = Template_helper_methods.scale(signal_area, signal1)
            f["ggH_0PM"] = (signal1, bins)
            
            self.signals[signal2_name] = signal2
            signal2, _ = np.histogram(signal2, bins=bins, range=(lowerlim, upperlim))
            signal2 = Template_helper_methods.scale(signal_area, signal2)
            f["ggH_0M"] = (signal2, bins)
            
            f["bkg_ggzz"] = self.scale_and_add_bkgs(bins, scaleTo=True)