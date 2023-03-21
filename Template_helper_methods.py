import os
import ROOT
import shutil
import uproot
import numpy as np
import argparse

def CombineParam(s):
    """Checks whether the argument you are passing is of form <Param Name>=<Param Numerical Value>
        Used as a type for argparse arguments i.e. parser.add_argument(--f, type=CombineParam)
        
    Parameters
    ----------
    s : str
        The argument being passed

    Returns
    -------
    str
        If the argument looks fine return the string as is

    Raises
    ------
    argparse.ArgumentTypeError
        Only one term should be on either side of the equals sign
    argparse.ArgumentTypeError
        Param value should be castable to a float
    """
    check = s.split('=')
    if len(check) != 2:
        raise argparse.ArgumentTypeError("Format should be <Param Name>=<Param Numerical Value>")
    
    try:
        float(check[1])
    except:
        raise argparse.ArgumentTypeError("Format should be <Param Name>=<Param Numerical Value>")
    
    return s


def scale(scaleto, counts, bins=[], return_scale_factor=False):
    """This function scales histograms according to their absolute area under the curve (no negatives allowed!)

    Parameters
    ----------
    scaleto : float
        The absolute area to scale to
    counts : list[Union[int,float]]
        A list of bin counts
    bins : list[float]
        The bins you want to use; use this option if you are passing a numpy histogram in (i.e. scale(1, *<numpy histogram>)), by default None
    return_scale_factor : bool
        Whether the scale factor being used will be returned
    
    Returns
    -------
    list[float]/Tuple(list[float], list[float])
        The scaled histogram bin counts or the scaled histogram, depending on whether you passed the bins in as well
    """
    counts = np.array(counts)
    counts = counts.astype(float)
    signs = np.sign(counts) #makes sure to preserve sign
    counts = np.abs(counts)
    
    if any(bins):
        if return_scale_factor:
            return signs*counts*scaleto/np.sum(counts), bins, scaleto/np.sum(counts)
        return signs*counts*scaleto/np.sum(counts), bins
    
    elif return_scale_factor:
        return signs*counts*scaleto/np.sum(counts), scaleto/np.sum(counts)
    
    return signs*counts*scaleto/np.sum(counts)

def Unroll_2D_OnShell(directory, fname):
    """Code written by Jeffrey Davis of happy hour cocktail fame to unroll a 2 dimensional histogram

    Parameters
    ----------
    directory : str
        The directory that you are inputting and outputting from
    fname : str
        The filename of what you are unrolling
    """
    if directory[-1] != '/':
        directory += '/'
    fname = fname.split('.')[0]
        
    histfile = ROOT.TFile.Open(directory+fname+'.root', "READ")
    fout = ROOT.TFile(directory+fname+'_unrolled.root',"RECREATE")
    fout.cd()
    
    for keyname in ['ggH_0PM', 'ggH_0M', 'bkg_ggzz']:
        hist = histfile.Get(keyname)
        
        xbins = hist.GetNbinsX()
        ybins = hist.GetNbinsY()

        temp_pos = ROOT.TH1F("temp_pos","",xbins*ybins,0,xbins*ybins)
        temp_neg = ROOT.TH1F("temp_neg","dif",xbins*ybins,0,xbins*ybins)
        #Unroll Hists
        indk = 0
        has_negative = False 
        for y in range (1,ybins+1):
            for x in range (1,xbins+1):
                binx_c = hist.GetXaxis().GetBinCenter(x)
                biny_c = hist.GetYaxis().GetBinCenter(y)
                ibin =  hist.FindBin(binx_c,biny_c)
                cont  = hist.GetBinContent(ibin)
                #put small values in empty background bins
                if cont == 0 : 
                    if "bkg" in hist.GetName():
                        intt = hist.Integral()
                        nb = ybins*xbins
                        contt = 0.1*intt*1.0/nb
                        # print ("found empty bin",contt)
                        hist.SetBinContent(ibin,contt)
                        # print (cont)
                if cont  < 0 :
                    has_negative = True
                    
        for y in range (1,ybins+1):
            for x in range (1,xbins+1):
                binx_c = hist.GetXaxis().GetBinCenter(x)
                biny_c = hist.GetYaxis().GetBinCenter(y)
                ibin =  hist.FindBin(binx_c,biny_c)
                cont  = hist.GetBinContent(ibin)
                if cont  < 0 :
                    temp_neg.Fill(indk,-1*cont)
                else :
                    temp_pos.Fill(indk,cont)
                    temp_pos.SetBinError(indk, np.sqrt(cont))
                indk = indk +1

        temp_name = hist.GetName()
        
        tpname = temp_name
        tnname = temp_name

        if (has_negative and ( "bkg" in tnname or "Data" in tnname  or "0PH" in tnname or "0PM" in tnname or "L1" in tnname or "0M" in tnname)):
            for y in range (1,ybins+1):
                for x in range (1,xbins+1):
                    binx_c = hist.GetXaxis().GetBinCenter(x)
                    biny_c = hist.GetYaxis().GetBinCenter(y)
                    ibin =  hist.FindBin(binx_c,biny_c)
                    cont  = hist.GetBinContent(ibin)

                #put small values in negtative background bins
                #Also put 0 in negative signal bins
                if cont  < 0 :
                    hist.SetBinContent(ibin,0)
                    print ("found negative bin",cont)
                    cont = 0
                if cont == 0 :
                    if "bkg" in hist.GetName():
                        intt = hist.Integral()
                        nb = ybins*xbins
                        contt = 0.1*intt*1.0/nb
                        print ("found empty bin",contt)
                        hist.SetBinContent(ibin,contt)
                        print (cont)
            
            temp_neg.SetName(tnname)
            temp_pos.SetName(tpname)

        elif (has_negative or not ( "bkg" in tnname or "Data" in tnname  or "0PH" in tnname or "0PM" in tnname or "L1" in tnname or "0M" in tnname) ):

            if "up" in tpname or "dn" in tpname :
                tpnm = tpname.split("_")
                tpnm.insert(2,"positive")
                tpname= tpnm[0]
                for ist in range(1,len(tpnm)):
                    tpname = tpname+"_"+tpnm[ist] 
            else :     
                tpname = tpname+"_positive"


            if "up" in tnname or "dn" in tnname :
                tnnm = tnname.split("_")
                tnnm.insert(2,"negative")
                tnname= tnnm[0]
                for ist in range(1,len(tnnm)):
                    tnname = tnname+"_"+tnnm[ist]  
            else :     
                tnname = tnname+"_negative"

                
            temp_neg.SetName(tnname)
            temp_pos.SetName(tpname)

        else:
        
            tnname = tnname.replace("0Xff_","0Mff_")
            tpname = tpname.replace("0Xff_","0Mff_")   
            
            temp_neg.SetName(tnname)
            temp_pos.SetName(tpname)
        if "data" in  tnname or "Data" in tnname : 
            
            temp_neg.SetName("data_obs")
            temp_pos.SetName("data_obs")
        
        
        
        if temp_pos.Integral() > 0:
            temp_pos.Write()
        if temp_neg.Integral() > 0:
            temp_neg.Write()
            
        print('Dumped Histogram into '+directory + fname+'_unrolled.root')