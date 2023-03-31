import numpy as np
import mplhep as hep
import matplotlib.pyplot as plt
plt.style.use(hep.style.ROOT)


def plot_overall_interference(terms, names,
                              output_directory, output_filename):
    """This function plots the overall plot of both interference and pure terms to plot everything

    Parameters
    ----------
    terms : list[tuple[Union[float, int]]]
        A list of all the pure sample and interference terms. This should be a list of (count, bin) pairs (i.e. numpy histograms)
    names : list[str]
        A list of the names for all of these terms
    output_directory : str
        The directory you would like to output to
    output_filename : str
        The filename you want to name the plots

    Returns
    -------
    tuple[numpy.ndarray, numpy.ndarray]
        the numpy histogram object of the overall sample
    """
    
    if output_directory[-1] != '/':
        output_directory += '/'
    
    plt.cla()
    plt.gca().axhline(zorder=-1, color='black', lw=2)
    
    bins = terms[0][1]
    overall = np.zeros(len(bins) - 1, dtype=float)
    for n, term in enumerate(terms):
        overall += term[0]
        if "_" in names[n]:
            hep.histplot(term, label=names[n], lw=2)
        else:
            hep.histplot(term, label=names[n], lw=4)
        
    # hep.histplot(overall, bins, label="Overall", lw=3, color='black')
    plt.legend(loc="upper right")
    plt.xlabel(r"$m_{4\mu}[GeV]$")
    plt.xlim(6,9)
    plt.tight_layout()
    plt.savefig(output_directory + output_filename + "overall_interference_plot.png")
    plt.savefig(output_directory + output_filename + "overall_interference_plot.pdf")
    
    return overall, bins