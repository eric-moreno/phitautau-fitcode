from hists_map import *
import numpy as np
import logging
from coffea import hist 
import hist as histhist

def intRegion(
    inhist: hist.Hist,
    theregion: str,
    nnCut: float,
    nnCut_loose: float,
    nnCut_fail: float,
    systematic: str = "nominal",
    samplelist: list = None,
    mslice: slice = None,
    mrebin: bool = False,
    debug: bool = False,
):
    """
    Return result of integrating over nnCuts.
    If there is a mass slice (mslice) also integrate over mass.
    """
    print("intRegion")
    print("theregion", theregion)
    print("nnCut", nnCut)
    print("nnCut_loose", nnCut_loose)
    print("nnCut_fail", nnCut_fail)
    print("systematic", systematic)
    print("samplelist", samplelist)
    # slice on nn
    if theregion == "pass":
        theslice = slice(nnCut, None)
        overflow_str = "over"
    elif theregion == "loosepass":
        theslice = slice(nnCut_loose, nnCut)
        overflow_str = "none"
    elif theregion == "fail":
        theslice = slice(nnCut_fail, nnCut_loose)
        overflow_str = "none"
    else:
        print("Unknown region", theregion)
        return
    the_int = inhist.integrate("nn_disc", theslice, overflow_str).integrate(
        "systematic", systematic
    )
    #print('THE INT', the_int)
    if True:
        print("nn slice ", theslice, "nn ", nnCut, " l ", nnCut_loose)
        print("overflow ", overflow_str, the_int.values())

    logging.debug(
        f"Values: {the_int.values()} with nncuts ({nnCut_loose} and {nnCut}) and overflow {overflow_str}"
    )

    # rebin in mass
    # if mrebin:
    #     original_bins = inhist.axis("massreg").edges()
    #     print(original_bins)
    #     print(inhist.axis("massreg").edges())
    #     # Define the new bin edges
    #     new_bin_edges = np.concatenate(([original_bins[0]], [original_bins[2]], [original_bins[4]], original_bins[6:]))
    #     print(new_bin_edges)
    #     # Perform rebinning
    #     the_int = the_int.rebin("massreg", hist.Bin("massreg", "massreg", new_bin_edges))
    #print(the_int)
        #the_int = the_int.rebin("massreg", hist.Bin("massreg", "massreg", mrebin))
    # if mrebin:
    #     the_int = custom_rebin(the_int, "massreg")

    # slice on mass
    if mslice is not None:
        if mslice.start is not None and mslice.stop is not None:
            overflow_str = "none"
        elif mslice.start is None and mslice.stop is not None:
            overflow_str = "under"
        elif mslice.start is not None and mslice.stop is None:
            overflow_str = "over"
        else:
            overflow_str = "allnan"
        the_int = the_int.integrate("massreg", mslice, overflow_str)
        if debug:
            print("massreg", overflow_str, the_int.values())
        logging.debug(
            f"Values: {the_int.values()} with mass slices and overflow {overflow_str}"
        )
    #print(the_int)
    # slice or sum samples
    if samplelist is not None:
        the_int = the_int.integrate("sample", samplelist).values(sumw2=True)[()]
    else:
        print("No samplelist provided, summing over all samples")
        the_int = the_int.sum("sample").values(sumw2=True)
        print("the_int", the_int)
        if () in the_int:
            the_int = the_int[()]
        else:
            mttbins = np.array(
            [
                0.0,
                10.0,
                20.0,
                30.0,
                40.0,
                50.0,
                60.0,
                70.0,
                80.0,
                90.0,
                100.0,
                110.0,
                120.0,
                130.0,
                140.0,
                150.0,
                200.0,
                250.0,
                300.0,
                350.0,
                400.0,
            ]
             )   
            #the_int = np.zeros(len(the_int), dtype=np.float32)
            the_int = (np.array([np.zeros(len(mttbins) - 1, dtype=np.float32)]))
    #print(the_int)
    if debug:
        print("debug", the_int)
    logging.debug(f"After integrating: {the_int}")
    print(the_int)
    return the_int


def getQCDFromData_old(
    inhist,
    region,
    nnCut,
    nnCut_loose,
    default=0.0,
    systematic="nominal",
    mslice=None,
    unc_scale=0.0,
    lowmassbin=2,
    highmassbin=-1,
    test = False
):
    """
    Get Data - other(MC) = QCD.
    :param default: sets event content for bins with 0 events
    :type default: float
    :param unc_scale: can force the errors to be at least some percent (potentially useful for bins with a lot of events)
    :type unc_scale: float
    """

    def clipZeros(data):
        data = np.clip(data, default, None)
        data[data == default] = 0.0
        return data

    qcd_data_full = intRegion(
        inhist["data_obs"],
        region,
        nnCut,
        nnCut_loose,
        systematic=systematic,
        mslice=mslice,
    )[0][lowmassbin:highmassbin]

    #print('qcd data full~!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    #print(qcd_data_full)

    other_mc = intRegion(
        inhist,
        region,
        nnCut,
        nnCut_loose,
        samplelist=[
            s.name
            for s in inhist.identifiers("sample")
            if s.name != "data_obs" and s.name != "multijet"
        ],
        systematic=systematic,
        mslice=mslice,
    )
    if test == True:
        print('inhist', inhist)
        print('region', region)
        print('samplelist', [
            s.name
            for s in inhist.identifiers("sample")
            if s.name != "data_obs" and s.name != "multijet"
        ])
        print('nnCut', nnCut)
        print('nnCut_loose', nnCut_loose)
        print('systematic', systematic)
        print('mslice', mslice)
    
    #print('other_mc~!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    #print(other_mc)

    #print('bins')
    #print(lowmassbin)
    #print(highmassbin)
    


    qcd_data = clipZeros(qcd_data_full - other_mc[0][lowmassbin:highmassbin])
    print("QCD DATA")
    print(qcd_data)
    # logging.debug(f"    Data: %s"%qcd_data_full)
    # logging.debug(f"    Other MC: %s"%other_mc[0][lowmassbin:highmassbin])
    logging.debug(f"    Data - MC: %s" % qcd_data)

    # QCD from data (data - other mc)
    qcd_temp = qcd_data

    # get variations ( from the poisson uncertainty on the number of events in each bin )
    qcd_data_altdn = clipZeros(
        qcd_data_full - other_mc[0][lowmassbin:highmassbin] * (1.0 + unc_scale)
    )
    qcd_data_altup = clipZeros(
        qcd_data_full - other_mc[0][lowmassbin:highmassbin] * (1.0 - unc_scale)
    )

    # QCD from data (data - other mc)
    qcd_data_w2 = np.clip(
        qcd_data_full - other_mc[1][lowmassbin:highmassbin], default, None
    )
    qcd_data_w2[qcd_data_w2 == default] = 0.0

    qcd_data_int = hist.poisson_interval(qcd_data, qcd_data_w2)
    #print('other mc')
    #print(other_mc)
    print('other mc 0')
    print(other_mc[0][lowmassbin:highmassbin])
    print('other mc 1')
    print(other_mc[1][lowmassbin:highmassbin])
    print('poisson interval')
    print(qcd_data_int)
    print('qcd_data')
    print(qcd_data)
    print('qcd_data_altdn')
    print(qcd_data_altdn)
    print('qcd_data_altup')
    print(qcd_data_altup)
    print('qcd_data_w2')
    print(qcd_data_w2)

    
    if test: 
        bins = np.array([20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 200, 250, 300, 350])
        import matplotlib.pyplot as plt
        bin_widths = np.diff(bins)
        bin_centers = bins[:-1] + bin_widths / 2

        # Plotting the histogram as a step plot
        plt.figure(figsize=(12, 8))
        plt.step(bin_centers, qcd_data_full, where='mid', label='data_full', linewidth=2)
        plt.step(bin_centers, other_mc[0][lowmassbin:highmassbin], where='mid', label='other_mc[0]', linewidth=2)
        plt.step(bin_centers, other_mc[1][lowmassbin:highmassbin], where='mid', label='other_mc[1]', linewidth=2)
        plt.step(bin_centers, qcd_data, where='mid', label='qcd data', linestyle='-.', linewidth=2)
        plt.step(bin_centers, qcd_data_w2, where='mid', label='qcd data w2', linestyle='-.', linewidth=2)
        plt.step(bin_centers, qcd_data_int[0], where='mid', label='poisson interval dn', linestyle='--', linewidth=2)
        plt.step(bin_centers, qcd_data_int[1], where='mid', label='poisson interval up', linestyle='--', linewidth=2)
        
        title = 'qcdcr_nom_fail'
        plt.xlabel('Bin Range')
        plt.ylabel('Counts')
        plt.title(title)
        plt.legend()
        plt.grid(True)
        plt.savefig(f'plots/{title}_hist.jpg')
    
    # qcd_temp_dn = np.minimum(
    #     np.array(
    #         [
    #             qcd_data_int[0][bi]
    #             if qcd_data_int[0][bi] >= 0.0 and not np.isnan(qcd_data_int[0][bi])
    #             else 0.0
    #             for bi in range(len(qcd_data))
    #         ]
    #     ),
    #     qcd_data_altdn,
    # )
    # qcd_temp_up = np.maximum(
    #     np.array(
    #         [
    #             qcd_data_int[1][bi]
    #             if qcd_data_int[1][bi] >= 0.0 and not np.isnan(qcd_data_int[1][bi])
    #             else 0.0
    #             for bi in range(len(qcd_data))
    #         ]
    #     ),
    #     qcd_data_altup,
    # )

    qcd_temp_dn = np.array([
        qcd_data_int[0][bi]
        if qcd_data_int[0][bi] >= 0.0 and not np.isnan(qcd_data_int[0][bi])
        else 0.0
        for bi in range(len(qcd_data))
        ])

    qcd_temp_up = np.array([
        qcd_data_int[1][bi]
        if qcd_data_int[1][bi] >= 0.0 and not np.isnan(qcd_data_int[1][bi])
        else 0.0
        for bi in range(len(qcd_data))
    ])


    return qcd_temp, qcd_temp_dn, qcd_temp_up

def getQCDFromData(
    inhist,
    region,
    nnCut,
    nnCut_loose,
    nnCut_fail,
    default=0.0,
    systematic="nominal",
    mslice=None,
    unc_scale=0.0,
    lowmassbin=2,
    highmassbin=-1,
    test=False,
    mrebin=False
):
    """
    Get Data - other(MC) = QCD.
    :param default: sets event content for bins with 0 events
    :type default: float
    :param unc_scale: can force the errors to be at least some percent (potentially useful for bins with a lot of events)
    :type unc_scale: float
    """

    def clipZeros(data):
        data = np.clip(data, default, None)
        data[data == default] = 0.0
        return data

    # Integrate data_obs histogram over the specified region and cuts
    data_obs_integrated = intRegion(
        inhist["data_obs"],
        region,
        nnCut,
        nnCut_loose,
        nnCut_fail,
        systematic=systematic,
        mslice=mslice,
        mrebin=mrebin
    )
    qcd_data_full = data_obs_integrated[0][lowmassbin:highmassbin]
    qcd_data_full_w2 = data_obs_integrated[1][lowmassbin:highmassbin]

    # Integrate all other MC samples except 'data_obs' and 'multijet'
    other_mc_integrated = intRegion(
        inhist,
        region,
        nnCut,
        nnCut_loose,
        nnCut_fail,
        samplelist=[
            s.name
            for s in inhist.identifiers("sample")
            if s.name != "data_obs" and s.name != "multijet"
        ],
        systematic=systematic,
        mslice=mslice,
        mrebin=mrebin
    )
    other_mc = other_mc_integrated[0][lowmassbin:highmassbin]
    other_mc_w2 = other_mc_integrated[1][lowmassbin:highmassbin]

    # Subtract other MC from data_obs and clip values to ensure no negative values
    qcd_data = clipZeros(qcd_data_full - other_mc)
    logging.debug(f"    Data - MC: %s" % qcd_data)

    # Calculate the Poisson interval for the data
    qcd_data_int = hist.poisson_interval(qcd_data_full, qcd_data_full_w2)
    other_mc_int = hist.poisson_interval(other_mc, other_mc_w2)

    # Calculate the up and down variations for data
    qcd_data_up = qcd_data_int[1]
    qcd_data_down = qcd_data_int[0]

    # Calculate the up and down variations for QCD
    qcd_temp_up = np.nan_to_num(clipZeros(qcd_data_up - other_mc_int[1]))
    qcd_temp_dn = np.nan_to_num(clipZeros(qcd_data_down - other_mc_int[0]))


    
    if 1: 
        bins = np.array([20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 200, 250, 300, 350])
        # print(bins.shape)
        # print(qcd_data_full.shape)
        # print(qcd_data.shape)
        # print(qcd_temp_dn.shape)
        # print(qcd_temp_up.shape)
        
        import matplotlib.pyplot as plt
        bin_widths = np.diff(bins)
        bin_centers = bins[:-1] + bin_widths / 2
        print(bin_centers.shape)
        # Plotting the histogram as a step plot
        plt.figure(figsize=(12, 8))
        plt.step(bin_centers, qcd_data_full, where='mid', label='data_full', linewidth=2)
        plt.step(bin_centers, other_mc, where='mid', label='other_mc', linewidth=2)
        plt.step(bin_centers, qcd_data, where='mid', label='qcd data', linestyle='-.', linewidth=2)
        plt.step(bin_centers, qcd_temp_up, where='mid', label='qcd temp up', linewidth=2)
        plt.step(bin_centers, qcd_temp_dn, where='mid', label='qcd temp dn', linewidth=2)
        
        title = test
        plt.xlabel('Bin Range')
        plt.ylabel('Counts')
        plt.title(title)
        plt.legend()
        plt.grid(True)
        plt.savefig(f'plots/{title}_hist.jpg')

    return qcd_data, qcd_temp_dn, qcd_temp_up


def getHist(
    h, var_name, lumifb, vars_cut, regionsel, blind, sigscale, rebin, debug=False
):
    """
    Get a histogram after integration and selections

    :param h: histogram
    :type: hist.Hist
    :param var_name: variable to keep
    :type str
    :param lumifb: luminosity
    :type float
    :param vars_cut:
    :param regionsel"
    :param blind:
    :param sigscale:
    :param rebin:

    """
    if debug:
        print(h)
    exceptions = ["process", var_name]
    exceptions.extend([var for var in vars_cut])
    if regionsel != "":
        exceptions.append("region")
    overflow_sum = "allnan"
    x = h.sum(
        *[ax for ax in h.axes() if ax.name not in exceptions], overflow=overflow_sum
    )
    for reg in regionsel:
        x = x.integrate("region", reg)

    for var, val in vars_cut.items():
        if len(val) == 0:
            continue
        if var != var_name:
            if debug:
                print("integrating ", var, val[0], val[1])
            if val[0] is not None and val[1] is not None:
                overflow_str = "none"
            elif val[0] is None and val[1] is not None:
                overflow_str = "under"
            elif val[0] is not None and val[1] is None:
                overflow_str = "over"
            else:
                overflow_str = "allnan"
            x = x.integrate(var, slice(val[0], val[1]), overflow=overflow_str)

    if var_name in vars_cut.keys():
        if debug:
            print(
                "integrating ", var_name, vars_cut[var_name][0], vars_cut[var_name][1]
            )
        x = x.rebin(
            var_name,
            hist.Bin(
                var_name,
                var_name,
                [
                    xe
                    for xe in x.axis(var_name).edges()
                    if (xe >= vars_cut[var_name][0] and xe <= vars_cut[var_name][1])
                ],
            ),
        )

    xaxis = var_name
    for ih, hkey in enumerate(x.identifiers("process")):
        x.identifiers("process")[ih].label = process_latex[hkey.name]

    if len(rebin) == 1:
        if rebin[0] > 1:
            print(rebin)
            x = x.rebin(xaxis, int(rebin[0]))
    else:
        # print("REBIN WTF")
        # print(rebin)
        # print('REBINNNNN')
        # sys.exit()
        x = x.rebin(xaxis, hist.Bin(xaxis, var_name, rebin))
    

    x_nobkg = x[nobkg]
    x_nosig = x[nosig]
    if debug:
        print("lumifb ", lumifb)
        print("sigscale ", sigscale)
    x_nosig.scale({p: lumifb for p in x_nosig.identifiers("process")}, axis="process")
    x_nobkg.scale(
        {p: lumifb * float(sigscale) for p in x_nobkg.identifiers("process")},
        axis="process",
    )
    x_data = x["data"]

    return x_nobkg + x_nosig + x_data


def coffea_to_hist(coffea_h, axis_name="massreg"):
    """
    Convert a 1D coffea histogram (with Weighted storage)
    into a hist.Hist object (boost-histogram-based).
    """
    # 1. Extract bin edges from the coffea axis
    edges = coffea_h.axis(axis_name).edges()

    # 2. Create an equivalent hist.Hist (with weight storage)
    h = histhist.Hist(
        hist.axis.Variable(edges),  # or Regular(...) if uniform
        storage=hist.storage.Weight()
    )

    # 3. Sum over any extra axes you don’t care about, e.g. "process"
    #    so that we're left with only the 1D axis.
    #    If you want to keep each process separate, you'd handle that differently.
    h_vals, h_vars = coffea_h.sum("process").values(sumw2=True)[()]  # 1D array

    # 4. Assign values and variances directly into the hist.Hist storage
    #    hist.Hist exposes a .view(flow=True/False) which gives you a buffer
    #    you can write to. For a 1D histogram, you can do:
    h.view(flow=False).value = h_vals
    h.view(flow=False).variance = h_vars

    return h

def hist_to_coffea(h_boost, axis_name="massreg"):
    """
    Convert a 1D hist.Hist (Weight storage)
    back into a coffea.hist.Hist.
    """

    # 1. Extract the edges from the hist.Hist
    edges = h_boost.axes[0].edges  # for a 1D hist, there's only one axis

    # 2. Create a new coffea histogram with these edges
    coffea_h = hist.Hist(
        "Events",
        hist.Bin(axis_name, axis_name, edges),
    )
    # By default, coffea uses arrays under the hood for ._sumw and ._sumw2

    # 3. Pull out values and variances from h_boost
    #    For 1D, h_boost.view(flow=False) is a 1D storage buffer
    vals = h_boost.view(flow=False).value
    vars = h_boost.view(flow=False).variance

    # 4. Assign them into coffea's internal arrays
    #    Typically, coffea uses some internal indexing, but for a plain 1D:
    for i in range(len(vals)):
        coffea_h._sumw[()][i] = vals[i]
        coffea_h._sumw2[()][i] = vars[i]

    return coffea_h
