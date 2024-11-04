import numpy as np
import logging
import rhalphalib as rl
from rhalphalib import MorphHistW2, AffineMorphTemplate
from utils import intRegion, getQCDFromData
from coffea import hist
import matplotlib.pyplot as plt

rl.ParametericSample.PreferRooParametricHist = True

debug_plots = False

def plot_histogram(data, bins, title):
    """
    Plots a histogram.
    
    Parameters:
    data (array): The counts in each bin.
    bins (array): The edges of the bins.
    title (str): The title for the histogram.
    """
    # Calculate the width of each bin for plotting
    bin_widths = np.diff(bins)
    bin_centers = bins[:-1] + bin_widths / 2

    # Plotting the histogram
    plt.figure(figsize=(10, 6))
    plt.bar(bin_centers, data, width=bin_widths, align='center', alpha=0.7, color='b')
    plt.xlabel('Bin Range')
    plt.ylabel('Counts')
    plt.title(title)
    plt.grid(True)
    plt.savefig(f'{title}_hist.jpg')

def plot_overlayed_histograms(nom_counts, shift_dn_counts, shift_up_counts, bins, title):
    """
    Plots overlayed histograms for nominal, shift down, and shift up data.

    Parameters:
    nom_counts (array): Counts for the nominal histogram.
    shift_dn_counts (array): Counts for the downward shifted histogram.
    shift_up_counts (array): Counts for the upward shifted histogram.
    bins (array): The edges of the bins.
    title (str): The title for the plot.
    """
    # Calculate the width of each bin for plotting
    bin_widths = np.diff(bins)
    bin_centers = bins[:-1] + bin_widths / 2

    # Plotting the histogram as a step plot
    plt.figure(figsize=(12, 8))
    plt.step(bin_centers, nom_counts, where='mid', label='Nominal', linewidth=2)
    plt.step(bin_centers, shift_dn_counts, where='mid', label='Shift Down', linestyle='--', linewidth=2)
    plt.step(bin_centers, shift_up_counts, where='mid', label='Shift Up', linestyle='-.', linewidth=2)

    plt.xlabel('Bin Range')
    plt.ylabel('Counts')
    plt.title(title)
    plt.legend()
    plt.grid(True)
    plt.savefig(f'plots/{title}_hist.jpg')

# Example bin edges and counts
bins = np.array([20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 200, 250, 300, 350])
#bins = np.array([20, 40, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 200, 250, 300, 350]) #wider bins lowmass problems


class Cards:
    def __init__(self, cat: str, year: str, no_syst: bool):
        self.logger = logging.getLogger("cards")

        self.__name = cat
        self.year = year
        self.islephad = cat != "hadhad"

        # regressed mass bins
        self.mttbins = np.array(
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

        self.mttbins_cut = np.array(
            [
                0.0,
                20.0,
                40.0,
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

        # bin slice
        self.lowbin_rebin = 1
        self.highbin_rebin = -1
        self.lowbin = 2
        self.highbin = -1
        self.mttbins_nom_rebin = self.mttbins_cut[self.lowbin_rebin : self.highbin_rebin]
        self.mttbins_nom = self.mttbins[self.lowbin : self.highbin]
        self.mttrange = range(len(self.mttbins_nom))

        # mass range
        self.mass = "massreg"
        self.massone = "massreg_one"
        self.lowqcdmass = 55.0 if self.islephad else 105.0
        self.highmass = 145.0 if self.islephad else 145.0
        self.lowqcdincrease = 0.5
        self.highbkgincrease = 0.3

        # pt range
        self.pt_min = 300

        # signal names
        self.doHtt = False
        if self.doHtt: 
            for i in range(100):
                print("DOING HTT DOING HTT DONT IGNORE ME")
        self.signame = "htt125" if self.doHtt else "phitt"
        self.sigexname = "htt125" if self.doHtt else "phitt50"
        #print(self.signame)
        
        # neural network cuts
        # pass region
        self.nnCut = {
            "hadhad": 0.9999,
            "hadel": 0.98,
            "hadmu": 0.98,
        }[cat]
        # loose-pass region
        self.nnCutLoose = {
            "hadhad": 0.995,
            "hadel": 0.9,
            "hadmu": 0.9,
        }[cat]

        self.nnCutFail = {
            "hadhad": None,
            "hadel": None,
            "hadmu": None,
        }[cat]
        
        # group samples
        self.sample_groups = {
            "data_obs": ["data"],
            "top": ["tt-had", "tt-semilep", "tt-dilep", "st"],
            "htt125": ["h125"],
            "multijet": ["qcd"],
            "dy": ["zem", "ztt"],
            "wlnu": ["wjets"],
            "vvqq": ["vv", "vqq"],
            "ignore": [],
        }

        # signal mass points
        self.masspoints = [
            "10",
            "20",
            "30",
            "40",
            "50",
            "75",
            "100",
            "125",
            "150",
            "200",
            "250",
            "300",
        ]
        for m in self.masspoints:
            if self.doHtt:
                self.sample_groups["ignore"].append(f"phi{m}")
            else:
                self.sample_groups[f"phitt{m}"] = [f"phi{m}"]


        #print(self.doHtt)
        #print(self.sample_groups["ignore"])
        #sys.exit()
        # met cuts
        cuts = {
            "lowmet": 50.0,
            "met": 75.0,
        }
        # met cut is tighter for hadhad region
        if not self.islephad:
            cuts["lowmet"] = 75.0
            cuts["met"] = 150.0

        slices = {
            "lowtohighmet": {
                "slice": slice(cuts["lowmet"], cuts["met"]),
                "overflow": "none",
            },
            "highmet": {"slice": slice(cuts["met"], None), "overflow": "over"},
            "lowmet": {"slice": slice(cuts["lowmet"], None), "overflow": "over"},
        }
        # exceptions for hadhad category
        if not self.islephad:
            slices["lowtohighmet"]["overflow"] = "under"

        # analysis regions
        self.regions = ["sig", "top_cr", "wlnu_cr", "qcd_cr"]

        # inverted regions (and their met cuts) for background estimation
        # nom: nominal region
        # fail: inverted region
        if self.islephad:
            # For LEPHAD: same as sig but with a low met cut
            self.invregions = {
                "nom": {
                    "sig": slices["highmet"],
                    "top_cr": slices["highmet"],
                    "wlnu_cr": slices["highmet"],
                    "qcd_cr": slices["highmet"],
                },
                "fail": {
                    "sig": slices["lowtohighmet"],
                    "top_cr": slices["lowtohighmet"],
                    "wlnu_cr": slices["lowtohighmet"],
                    "qcd_cr": slices["lowtohighmet"],
                },
            }
        else:
            # For HADHAD: named "faildphi" or "dphi_inv" in the histogram names
            # https://github.com/drankincms/boostedhiggs/blob/c21f38fcec3b0df75a4327f2fd04139a840adbd5/boostedhiggs/httprocessor.py#L649
            self.invregions = {
                "nom": {
                    "sig": slices["highmet"],
                    "top_cr": slices["highmet"],
                    "wlnu_cr": slices["lowmet"],
                    "qcd_cr": slices["lowmet"],
                },
                "fail": {
                    "sig": slices["highmet"],
                    "top_cr": slices["highmet"],
                    "wlnu_cr": slices["lowmet"],
                    "qcd_cr": slices["lowmet"],
                },
            }

        # nn regions
        self.nnregions = ["fail", "loosepass", "pass"]

        # variation of background yields for systematics
        self.variations = ["nom", "dn", "up"]#, "dn2", "up2"]

        # initialize model
        self.model = rl.Model(f"{cat}Model")

        # initialize observable
        self.mtt = rl.Observable(self.mass, self.mttbins_nom)
        self.mttrebin = rl.Observable(self.mass, self.mttbins_nom_rebin)
        self.mttone = rl.Observable(
            self.massone,
            np.array([self.mttbins[self.lowbin], self.mttbins[self.highbin - 1]]),
        )

        # initialize nuisances
        self.nuisances(no_syst)

    def nuisances(self, no_syst: bool):
        self.logger.info(f"Initializing nuisances")
        cat = self.__name

        #add UE and JES/JER
        # syst_dict = {
        #     samp: {
        #         "uescale": [
        #             rl.NuisanceParameter("uescale", "shape"),
        #             "UESDown",
        #             "UESUp",
        #         ],
        #         "jescale": [
        #             rl.NuisanceParameter("jescale", "shape"),
        #             "JESDown",
        #             "JESUp",
        #         ],
        #         "jeresol": [
        #             rl.NuisanceParameter("jeresol", "shape"),
        #             "JERDown",
        #             "JERUp",
        #         ],
        #     }
        #     for samp in self.sample_groups
        #     if samp not in ["data_obs", "ignore", "multijet"]
        # }

        syst_dict = { 
                samp: {}
            for samp in self.sample_groups 
            if samp not in ["data_obs", "ignore", "multijet"]
        }

        lumi_list = {
            "2016APV": {"lumi_16APV": rl.NuisanceParameter("CMS_lumi_16APV", "lnN"), "lumi_all": rl.NuisanceParameter("CMS_lumi_all", "lnN")},
            "2016": {"lumi_16": rl.NuisanceParameter("CMS_lumi_16", "lnN"), "lumi_all": rl.NuisanceParameter("CMS_lumi_all", "lnN")},
            "2017": {"lumi_17": rl.NuisanceParameter("CMS_lumi_17", "lnN"), "lumi_all": rl.NuisanceParameter("CMS_lumi_all", "lnN"), "lumi_1718": rl.NuisanceParameter("CMS_lumi_1718", "lnN")},
            "2018": {"lumi_18": rl.NuisanceParameter("CMS_lumi_18", "lnN"), "lumi_all": rl.NuisanceParameter("CMS_lumi_all", "lnN"), "lumi_1718": rl.NuisanceParameter("CMS_lumi_1718", "lnN")},
        }[self.year]
        
        syst_dict_JER = {    
            "uescale": [
                rl.NuisanceParameter("uescale", "lnN"),
                "UESDown",
                "UESUp",
            ],
            "jescale": [
                rl.NuisanceParameter("jescale", "lnN"),
                "JESDown",
                "JESUp",
            ],
            "jeresol": [
                rl.NuisanceParameter("jeresol", "lnN"),
                "JERDown",
                "JERUp",
            ],
        }

        self.syst_dict_JER = syst_dict_JER

        # syst_dict_JER = {    
        #     "uescale":
        #         rl.NuisanceParameter("uescale", "shape"),
        #     "jescale":
        #         rl.NuisanceParameter("jescale", "shape"),
        #     "jeresol":
        #         rl.NuisanceParameter("jeresol", "shape"),
        #     "l1prefire": 
        #         rl.NuisanceParameter("l1prefire", "shape"),
        # }

        syst_dict_norm = {
            "CMS_vvqq_norm": rl.NuisanceParameter("CMS_vvqq_norm", "lnN"),
            "CMS_top_norm": rl.NuisanceParameter("CMS_top_norm", "lnN"),
            "CMS_wlnu_norm": rl.NuisanceParameter("CMS_wlnu_norm", "lnN") ,
            "CMS_dy_norm": rl.NuisanceParameter("CMS_dy_norm", "lnN"),
            f"CMS_trig_{cat}": rl.NuisanceParameter(f"CMS_trig_{cat}", "lnN"),
            f"CMS_id_{cat}": rl.NuisanceParameter(f"CMS_id_{cat}", "lnN"),
        }
        
        if self.year not in ["2018"]:
            syst_dict_prefire = {"l1prefire": [
                rl.NuisanceParameter("l1prefire", "lnN"),
                "L1PreFiringDown",
                "L1PreFiringUp",
            ]}
        else:
            syst_dict_prefire = {}
        #syst_dict_prefire = {}


        syst_dict_top = {"toppt": [
            rl.NuisanceParameter("toppt", "lnN"),
           "nominal",
           "TopPtReweightUp",]
        }
        # add prefire uncertainty
        # if self.year not in ["2018"]:
        #     for sample in syst_dict:
        #         syst_dict[sample]["l1prefire"] = [
        #             rl.NuisanceParameter("l1prefire", "shape"),
        #             "L1PreFiringDown",
        #             "L1PreFiringUp",
        #         ]
        
        syst_dict_cat = syst_dict 
        syst_dict_UPDOWN = syst_dict_JER | syst_dict_prefire 
        syst_dict = syst_dict_norm | lumi_list | syst_dict_top
        # syst_dict["top"] = [
        #    rl.NuisanceParameter("toppt", "shape"),
        #    "nominal",
        #    "TopPtReweightUp",
        # ]
        # add toppt uncertainty for top region
        
        # syst_dict["top"]["toppt"] = [
        #       rl.NuisanceParameter("toppt", "shape"),
        #       "nominal",
        #       "TopPtReweightUp",
        # ]

        #print(syst_dict)
        #print(syst_dict.keys())
        #sys.exit()
        # dictionary of systematic uncertainty variations to build from existing mass templates
        if no_syst:
            self.syst_dict = {}
            self.syst_dict_cat = syst_dict_cat
            self.syst_dict_UPDOWN = {}
        else:
            self.syst_dict = syst_dict
            self.syst_dict_cat = syst_dict_cat
            self.syst_dict_UPDOWN = syst_dict_UPDOWN

        # QCD normalization
        self.qcdnormSF = rl.IndependentParameter(f"qcdnormSF_{cat}", 1.0, 0, 10)
        self.qcdnormSF_top = rl.IndependentParameter(f"qcdnormSF_top_{cat}", 1.0, 0, 10)
        self.qcdnormSF_wlnu = rl.IndependentParameter(
            f"qcdnormSF_wlnu_{cat}", 1.0, 0, 10
        )

        # Bkg efficiency
        self.bkgeffSF = rl.IndependentParameter(f"bkgeffSF_{cat}", 1.0, 0, 10)
        self.bkgLeffSF = rl.IndependentParameter(f"bkgLeffSF_{cat}", 1.0, 0, 10)

        # Top efficiency
        self.topeffSF = rl.IndependentParameter(f"topeffSF_{cat}", 1.0, 0, 10)
        self.topLeffSF = rl.IndependentParameter(f"topLeffSF_{cat}", 1.0, 0, 10)

        # WJets efficiency
        self.wlnueffSF = rl.IndependentParameter(f"wlnueffSF_{cat}", 1.0, 0, 10)
        self.wlnuLeffSF = rl.IndependentParameter(f"wlnuLeffSF_{cat}", 1.0, 0, 10)

        self.rh125 = rl.NuisanceParameter(f"r_h125_{cat}", "lnN")

        # DY+Jets efficiency
        self.dy_eff = rl.IndependentParameter(f"dy_eff_{cat}", 1.0, 0, 10)

        # top and wlnu shape
        self.top_highmass = [
            rl.NuisanceParameter(f"top_highmass_bin{ix}_{cat}", "shape")
            for ix in self.mttrange
            if self.mttbins_nom[ix] > self.highmass
        ]
        self.wlnu_highmass = [
             rl.NuisanceParameter(f"wlnu_highmass_bin{ix}_{cat}", "shape")
             for ix in self.mttrange
             if self.mttbins_nom[ix] > self.highmass
         ]
        highmassx = -1
        for ix in self.mttrange:
            if self.mttbins_nom[ix] > self.highmass:
                if highmassx == -1:
                    highmassx = ix
        self.highmassx = highmassx
        
        # QCD shape
        self.qcd_fail_1 = rl.NuisanceParameter(f"qcd_Rfail_method1_{cat}", "shape")
        self.qcd_loosepass_1 = rl.NuisanceParameter(f"qcd_Rloosepass_method1_{cat}", "shape")
        self.qcd_pass_1 = rl.NuisanceParameter(f"qcd_Rpass_method1_{cat}", "shape")

        self.qcd_fail_2 = rl.NuisanceParameter(f"qcd_Rfail_method2_{cat}", "shape")
        self.qcd_loosepass_2 = rl.NuisanceParameter(f"qcd_Rloosepass_method2_{cat}", "shape")
        self.qcd_pass_2 = rl.NuisanceParameter(f"qcd_Rpass_method2_{cat}", "shape")

        # self.qcd_fail_1 = rl.NuisanceParameter(f"qcd_Rfail_method1_{cat}", "lnN")
        # self.qcd_loosepass_1 = rl.NuisanceParameter(f"qcd_Rloosepass_method1_{cat}", "lnN")
        # self.qcd_pass_1 = rl.NuisanceParameter(f"qcd_Rpass_method1_{cat}", "lnN")

        # self.qcd_fail_2 = rl.NuisanceParameter(f"qcd_Rfail_method2_{cat}", "lnN")
        # self.qcd_loosepass_2 = rl.NuisanceParameter(f"qcd_Rloosepass_method2_{cat}", "lnN")
        # self.qcd_pass_2 = rl.NuisanceParameter(f"qcd_Rpass_method2_{cat}", "lnN")

        # mass scale 
        #self.m_scale = rl.NuisanceParameter(f"massscale_{cat}", "shape")
        #self.m_scale_bkg = rl.NuisanceParameter(f"massscale_bkg_{cat}", "shape")

        self.lowmass_combined = False
        # add uncertainty for low mass

        if not self.lowmass_combined: 
            self.qcd_lowmass_1 = [
                rl.NuisanceParameter(f"qcd_lowmass_method1_bin{ix}_{cat}", "shape")
                for ix in self.mttrange
                if self.mttbins_nom[ix] < self.lowqcdmass
            ]

            self.qcd_lowmass_2 = [
                rl.NuisanceParameter(f"qcd_lowmass_method2_bin{ix}_{cat}", "shape")
                for ix in self.mttrange
                if self.mttbins_nom[ix] < self.lowqcdmass
            ]

            self.qcd_lowmass_top_1 = [
                rl.NuisanceParameter(f"qcd_lowmass_top_method1_bin{ix}_{cat}", "shape")
                for ix in self.mttrange
                if self.mttbins_nom[ix] < self.lowqcdmass
            ]
            self.qcd_lowmass_top_2 = [
                rl.NuisanceParameter(f"qcd_lowmass_top_method2_bin{ix}_{cat}", "shape")
                for ix in self.mttrange
                if self.mttbins_nom[ix] < self.lowqcdmass
            ]

            self.qcd_lowmass_wlnu_1 = [
                rl.NuisanceParameter(f"qcd_lowmass_wlnu_method1_bin{ix}_{cat}", "shape")
                for ix in self.mttrange
                if self.mttbins_nom[ix] < self.lowqcdmass
            ]
            self.qcd_lowmass_wlnu_2 = [
                rl.NuisanceParameter(f"qcd_lowmass_wlnu_method2_bin{ix}_{cat}", "shape")
                for ix in self.mttrange
                if self.mttbins_nom[ix] < self.lowqcdmass
            ]

        elif self.lowmass_combined:
            #combined lowmass
            self.qcd_lowmass_1 = [
                rl.NuisanceParameter(f"qcd_lowmass_method1_comb", "shape")
            ]

            self.qcd_lowmass_2 = [
                rl.NuisanceParameter(f"qcd_lowmass_method2_comb", "shape")
            ]

            self.qcd_lowmass_top_1 = [
                rl.NuisanceParameter(f"qcd_lowmass_top_method1_comb", "shape")
            ]
            self.qcd_lowmass_top_2 = [
                rl.NuisanceParameter(f"qcd_lowmass_top_method2_comb", "shape")
            ]

            self.qcd_lowmass_wlnu_1 = [
                rl.NuisanceParameter(f"qcd_lowmass_wlnu_method1_comb", "shape")
            ]
            self.qcd_lowmass_wlnu_2 = [
                rl.NuisanceParameter(f"qcd_lowmass_wlnu_method2_comb", "shape")
        ]

        self.sys_smear = rl.NuisanceParameter('CMS_resonance_smear', 'shape')
        self.sys_shift = rl.NuisanceParameter('CMS_resonance_shift', 'shape')
        self.sys_smear_nonres = rl.NuisanceParameter('CMS_nonresonance_smear', 'shape')
        self.sys_shift_nonres = rl.NuisanceParameter('CMS_nonresonance_shift', 'shape')
        # self.qcd_lowmass = rl.NuisanceParameter(f"qcd_lowmass_bin_{cat}", "shape")
        # self.qcd_lowmass_top = rl.NuisanceParameter(f"qcd_lowmass_top_bin_{cat}", "shape")
        # self.qcd_lowmass_wlnu = rl.NuisanceParameter(f"qcd_lowmass_wlnu_bin_{cat}", "shape")
        
    def get_qcd(self, hists: dict, region: str, default: float = 0.0, test=False):
        """
        Get QCD estimate by subtracting other_MC from data.

        :param hists: dictionary with "data_obs" and "otherMC" histograms
        :type hists: dict
        :param region: string: pass/loosepass/fail
        :type region: str
        :param default: sets event content for bins with 0 events
        :type default: float
        """
        return getQCDFromData(
            hists, region, self.nnCut, self.nnCutLoose, self.nnCutFail, default=default, test=test, mrebin=False
        )

    def _get_region(
        self, h: hist.Hist, region: str, syst: str = "nominal", debug: bool = False, rebin = False
    ):
        """
        Get region by integrating over nn cut and/or met/mass/pt cuts
        """
        return intRegion(
            h,
            region,
            self.nnCut,
            self.nnCutLoose,
            self.nnCutFail,
            systematic=syst,
            debug=debug,
            mrebin = rebin
        )

    def _events(self, template, singlebin=False, clip=True):
        """
        Get numpy array of event yields
        """
        templ = template[0][self.lowbin : self.highbin]
        if singlebin:
            templ = np.array([np.sum(templ)])
        if clip:
            templ = np.clip(templ, 0.0, None)
        # remove negative events
        for iv, val in enumerate(templ):
            if val < 0.0:
                templ[iv] = 0.0
        return templ

    def _sumw(self, template, singlebin=False):
        """
        Get sum of weights
        """
        templ = template[1][self.lowbin : self.highbin]
        if singlebin:
            templ = np.array([np.sum(templ)])
        return templ

    def get_template(self, h, region, syst, singlebin, rebin, debug, clip=True):
        """
        Get template ntuple
        """
        tempint = self._get_region(h, region, syst, debug, rebin=rebin)
        events = self._events(tempint, singlebin, clip)
        sumw = self._sumw(tempint, singlebin)
        template = (
            events,
            (self.mttone.binning if singlebin else (self.mttrebin.binning if rebin else self.mtt.binning)),
            (self.mttone.name if singlebin else (self.mttrebin.binning if rebin else self.mtt.name)),
            sumw,
        )
        return template, events
    

    def systs_mass(self, h, region, sample, sample_template, template_nonrl, singlebin=False, analysis_region=None,):
        """
        Add mass systematics
        """
        cat = self.__name
        syst_template = sample_template
        
        if (
            sample.name == "htt125" or sample.name == "dy" or "phi" in sample.name
        ) and not singlebin:
            # why not use rl mass shift functions?

            nominal = self._get_region(h, region)
            nom = self._events(nominal, clip=False)
            nom_full = nominal[0]
            shift_dn = nom_full[self.lowbin - 1 : self.highbin - 1]
            shift_up = nom_full[
                self.lowbin + 1 : self.highbin + 1 if self.highbin != -1 else None
            ]
            shiftwidth = self.mttbins[1] - self.mttbins[0]

            bins = np.array([20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 200, 250, 300, 350])

            if debug_plots:
                plot_histogram(nom, bins, f'{sample.name}_NominalHistogram')
                plot_histogram(shift_dn, bins, f'{sample.name}_ShiftDownHistogram')
                plot_histogram(shift_up, bins, f'{sample.name}_ShiftUpHistogram')
                plot_overlayed_histograms(nom, shift_dn, shift_up, bins, F'{sample.name}_OverlayedHistograms')

            # dnfrac = np.array(
            #     [
            #         shiftwidth
            #         / (
            #             self.mttbins[self.lowbin + ib]
            #             - self.mttbins[self.lowbin + ib - 1]
            #         )
            #         for ib in range(len(nom))
            #     ]
            # )
            # shift_dn = shift_dn * dnfrac + nom * (1.0 - dnfrac)

            # upfrac = np.array(
            #     [
            #         shiftwidth
            #         / (
            #             self.mttbins[self.lowbin + ib + 2]
            #             - self.mttbins[self.lowbin + ib + 1]
            #         )
            #         for ib in range(len(nom))
            #     ]
            # )
            # shift_up = shift_up * upfrac + nom * (1.0 - upfrac)
            shiftfrac = np.array([shiftwidth/(self.mttbins[self.lowbin+ib+1]-self.mttbins[self.lowbin+ib]) for ib in range(len(nom_full[self.lowbin:self.highbin]))]) # this accounts for variable bin widths
            shiftfrac_dn = np.insert(shiftfrac, 0, shiftfrac[0])[:-1]
            shiftfrac_up = np.append(shiftfrac, shiftfrac[-1])[1:]
            shift_dn = shift_dn*shiftfrac_dn + nom_full[self.lowbin:self.highbin]*(1.-shiftfrac_dn)
            shift_up = shift_up*shiftfrac_up + nom_full[self.lowbin:self.highbin]*(1.-shiftfrac_up)
            if debug_plots: 
                plot_histogram(shift_dn, bins, f'{sample.name}_ShiftDownHistogramAfterMultiplication')
                plot_histogram(shift_up, bins, f'{sample.name}_ShiftUpHistogramAfterMultiplication')
                plot_overlayed_histograms(nom, shift_dn, shift_up, bins, f'{analysis_region}_{region}_{sample.name}_MassShift')

            # syst_template.setParamEffect(
            #     self.m_scale,
            #     np.divide(shift_dn, nom_full[self.lowbin:self.highbin], out=np.ones_like(nom_full[self.lowbin:self.highbin]), where=nom_full[self.lowbin:self.highbin] > 0.0),
            #     np.divide(shift_up, nom_full[self.lowbin:self.highbin], out=np.ones_like(nom_full[self.lowbin:self.highbin]), where=nom_full[self.lowbin:self.highbin] > 0.0),
            # )

            # print(np.divide(shift_dn, nom_full[self.lowbin:self.highbin], out=np.ones_like(nom_full[self.lowbin:self.highbin]), where=nom_full[self.lowbin:self.highbin] > 0.0))
            # print(np.divide(shift_up, nom_full[self.lowbin:self.highbin], out=np.ones_like(nom_full[self.lowbin:self.highbin]), where=nom_full[self.lowbin:self.highbin] > 0.0))
            # sys.exit()
        if (
            sample.name == "top" or sample.name == "wlnu"
        ) and not singlebin:
            # why not use rl mass shift functions?

            nominal = self._get_region(h, region)
            nom = self._events(nominal, clip=False)
            nom_full = nominal[0]
            shift_dn = nom_full[self.lowbin - 1 : self.highbin - 1]
            shift_up = nom_full[
                self.lowbin + 1 : self.highbin + 1 if self.highbin != -1 else None
            ]
            shiftwidth = self.mttbins[1] - self.mttbins[0]


            bins = np.array([20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 200, 250, 300, 350])

            if debug_plots: 
                plot_overlayed_histograms(nom, shift_dn, shift_up, bins, F'{cat}_{region}_{sample.name}_OverlayedHistograms')

            # dnfrac = np.array(
            #     [
            #         shiftwidth
            #         / (
            #             self.mttbins[self.lowbin + ib]
            #             - self.mttbins[self.lowbin + ib - 1]
            #         )
            #         for ib in range(len(nom))
            #     ]
            # )
            # shift_dn = shift_dn * dnfrac + nom * (1.0 - dnfrac)

            # upfrac = np.array(
            #     [
            #         shiftwidth
            #         / (
            #             self.mttbins[self.lowbin + ib + 2]
            #             - self.mttbins[self.lowbin + ib + 1]
            #         )
            #         for ib in range(len(nom))
            #     ]
            # )
            # shift_up = shift_up * upfrac + nom * (1.0 - upfrac)
            shiftfrac = np.array([shiftwidth/(self.mttbins[self.lowbin+ib+1]-self.mttbins[self.lowbin+ib]) for ib in range(len(nom_full[self.lowbin:self.highbin]))]) # this accounts for variable bin widths
            shiftfrac_dn = np.insert(shiftfrac, 0, shiftfrac[0])[:-1]
            shiftfrac_up = np.append(shiftfrac, shiftfrac[-1])[1:]
            shift_dn = shift_dn*shiftfrac_dn + nom_full[self.lowbin:self.highbin]*(1.-shiftfrac_dn)
            shift_up = shift_up*shiftfrac_up + nom_full[self.lowbin:self.highbin]*(1.-shiftfrac_up)

            if debug_plots:
                plot_overlayed_histograms(nom, shift_dn, shift_up, bins, f'{analysis_region}_{region}_{sample.name}_MassShift')

            # syst_template.setParamEffect(
            #     self.m_scale_bkg,
            #     np.divide(shift_dn, nom_full[self.lowbin:self.highbin], out=np.ones_like(nom_full[self.lowbin:self.highbin]), where=nom_full[self.lowbin:self.highbin] > 0.0),
            #     np.divide(shift_up, nom_full[self.lowbin:self.highbin], out=np.ones_like(nom_full[self.lowbin:self.highbin]), where=nom_full[self.lowbin:self.highbin] > 0.0),
            # )
        

        smear_syst=True

        if smear_syst and sample.name not in ['top', 'multijet', 'wlnu', 'data_obs']:
            mtempl = MorphHistW2(template_nonrl)
            _up = mtempl.get(smear=1 + 0.20)
            _down = mtempl.get(smear=1 - 0.20)
            syst_template.setParamEffect(self.sys_smear, _up[:-1], _down[:-1], scale=1)

            _up = mtempl.get(shift = 2)
            _down = mtempl.get(shift = -2)
            syst_template.setParamEffect(self.sys_shift, _up[:-1], _down[:-1], scale=1)
        
        if smear_syst and sample.name in ['top', 'wlnu']:
            mtempl = MorphHistW2(template_nonrl)
            _up = mtempl.get(smear=1 + 0.20)
            _down = mtempl.get(smear=1 - 0.20)
            syst_template.setParamEffect(self.sys_smear_nonres, _up[:-1], _down[:-1], scale=1)

            _up = mtempl.get(shift = 2)
            _down = mtempl.get(shift = -2)
            syst_template.setParamEffect(self.sys_shift_nonres, _up[:-1], _down[:-1], scale=1)

        return syst_template

    def systs_shape(self, h, region, sample, sample_template, singlebin=False):
        """
        Add systematics with shape samples
        """

        syst_template = sample_template
        #print(region)
        #print(sample.name)
        #print(self.syst_dict)
        #print('done')
        if sample.name in ["data_obs", "ignore", "multijet"]: 
            pass
        else:
            tempint = self._get_region(h, region)
            nom = self._events(tempint)

            for syst in self.syst_dict_UPDOWN:
                nuisance, syst_dn, syst_up = self.syst_dict_UPDOWN[syst]
                up = self._events(self._get_region(h, region, syst_up), clip=False)
                dn = self._events(self._get_region(h, region, syst_dn))
                up_var = (
                    np.divide(up, nom, out=np.ones_like(nom), where=nom > 0.0)
                    if (up != nom).all()
                    else np.ones_like(nom) * 1.001
                )
                dn_var = (
                    np.divide(dn, nom, out=np.ones_like(nom), where=nom > 0.0)
                    if (dn != nom).all()
                    else np.ones_like(nom) * 0.999
                )
                syst_template.setParamEffect(
                    nuisance,
                    np.array([np.sum(up_var)]) if singlebin else up_var,
                    np.array([np.sum(dn_var)]) if singlebin else dn_var,
                )

            if sample.name == "top":
                nuisance, syst_dn, syst_up = rl.NuisanceParameter("toppt", "shape"), "nominal", "TopPtReweightUp"
                #nuisance, syst_dn, syst_up = self.syst_dict["toppt"]

                up = self._events(self._get_region(h, region, syst_up), clip=False)
                dn = self._events(self._get_region(h, region, syst_dn))

                up_var = (
                    np.divide(up, nom, out=np.ones_like(nom), where=nom > 0.0)
                    if (up != nom).all()
                    else np.ones_like(nom) * 1.001
                )
                dn_var = (
                    np.divide(dn, nom, out=np.ones_like(nom), where=nom > 0.0)
                    if (dn != nom).all()
                    else np.ones_like(nom) * 0.999
                )
            
                # syst_template.setParamEffect(
                #     nuisance,
                #     np.array([np.sum(up_var)]) if singlebin else up_var,
                #     np.array([np.sum(dn_var)]) if singlebin else dn_var,
                # )
        #print(syst_template)
        #print('syst template')
        return syst_template

    def systs_norm(self, region, events, qcd, sample, sample_template, singlebin=False, analysis_region=None):
        """
        Add systematics for normalization
        """

        cat = self.__name
        syst_template = sample_template

        if sample.name == "top":
            syst_template.setParamEffect(
                self.syst_dict["CMS_top_norm"], 1.10
            )
            if not singlebin:
                for imx in range(len(self.top_highmass)):
                    syst_template.setParamEffect(
                        self.top_highmass[imx],
                        np.array(
                            [
                                1.0 - self.highbkgincrease
                                if ix == imx + self.highmassx
                                else 1.0
                                for ix in range(len(events))
                            ]
                        ),
                        np.array(
                            [
                                1.0 + self.highbkgincrease
                                if ix == imx + self.highmassx
                                else 1.0
                                for ix in range(len(events))
                            ]
                        ),
                    )
                
        if sample.name == "wlnu":
            syst_template.setParamEffect(
                self.syst_dict["CMS_wlnu_norm"], 1.10
            )
            if not singlebin:
                for imx in range(len(self.wlnu_highmass)):
                    syst_template.setParamEffect(
                        self.wlnu_highmass[imx],
                        np.array(
                            [
                                1.0 - self.highbkgincrease
                                if ix == imx + self.highmassx
                                else 1.0
                                for ix in range(len(events))
                            ]
                        ),
                        np.array(
                            [
                                1.0 + self.highbkgincrease
                                if ix == imx + self.highmassx
                                else 1.0
                                for ix in range(len(events))
                            ]
                        ),
                    )

        # if sample.name == "vvqq":
        #     syst_template.setParamEffect(
        #         rl.NuisanceParameter("CMS_vvqq_norm", "lnN"), 1.20
        #     )
        # if sample.name == "dy":
        #     syst_template.setParamEffect(
        #         rl.NuisanceParameter("CMS_vvqq_norm", "lnN"), 1.05
        #     )

        if sample.name == "vvqq":
            syst_template.setParamEffect(
                self.syst_dict["CMS_vvqq_norm"], 1.20
            )
        if sample.name == "dy":
            syst_template.setParamEffect(
                self.syst_dict["CMS_dy_norm"], 1.05
            )


        # lumi_16 = rl.NuisanceParameter("CMS_lumi_16", "lnN")
        # lumi_16APV = rl.NuisanceParameter("CMS_lumi_16APV", "lnN")
        # lumi_17 = rl.NuisanceParameter("CMS_lumi_17", "lnN")
        # lumi_18 = rl.NuisanceParameter("CMS_lumi_18", "lnN")
        # lumi_all = rl.NuisanceParameter("CMS_lumi_all", "lnN")
        # lumi_1718 = rl.NuisanceParameter("CMS_lumi_1718", "lnN")
        lumi_list = {
            "2016APV": ["lumi_16APV", "lumi_all"],
            "2016": ["lumi_16", "lumi_all"],
            "2017": ["lumi_17", "lumi_all", "lumi_1718"],
            "2018": ["lumi_18", "lumi_all", "lumi_1718"],
        }[self.year]
        lumi_vals = {
            "2016APV": [1.01, 1.006],
            "2016": [1.01, 1.006],
            "2017": [1.02, 1.009, 1.006],
            "2018": [1.015, 1.02, 1.002],
        }[self.year]
        
        if sample.name != "multijet":
            syst_template.setParamEffect(
                self.syst_dict[f"CMS_trig_{cat}"], 1.02
            )
            if self.islephad:
                syst_template.setParamEffect(
                    self.syst_dict[f"CMS_id_{cat}"], 1.02
                )
            for il, lumi in enumerate(lumi_list):
                syst_template.setParamEffect(self.syst_dict[lumi], lumi_vals[il])

            
                    
        else:
            if not singlebin:
                qcd_shape_dn_1 = np.divide(
                    qcd["dn"],
                    qcd["nom"],
                    out=np.ones_like(qcd["nom"]),
                    where=qcd["nom"] > 0.0,
                )
                qcd_shape_up_1 = np.divide(
                    qcd["up"],
                    qcd["nom"],
                    out=np.ones_like(qcd["nom"]),
                    where=qcd["nom"] > 0.0,
                )

                if debug_plots: 
                    plot_overlayed_histograms(qcd["nom"], qcd["dn"],  qcd["up"], bins, f'{analysis_region}_{region}_{sample.name}_Method1')
                    plot_overlayed_histograms(qcd["nom"], qcd["dn2"],  qcd["up2"], bins, f'{analysis_region}_{region}_{sample.name}_Method2')

                if "dn2" in qcd.keys():
                    qcd_shape_dn_2 = np.divide(
                        qcd["dn2"],
                        qcd["nom"],
                        out=np.ones_like(qcd["nom"]),
                        where=qcd["nom"] > 0.0,
                    )
                if "up2" in qcd.keys():
                    qcd_shape_up_2 = np.divide(
                        qcd["up2"],
                        qcd["nom"],
                        out=np.ones_like(qcd["nom"]),
                        where=qcd["nom"] > 0.0,
                    )
                
                def symmetrize_binwise(dn, up):
                    sym_dn = np.copy(dn)
                    sym_up = np.copy(up)
                    for i in range(len(dn)):
                        fluct_dn = 1.0 - dn[i]
                        fluct_up = up[i] - 1.0
                        min_fluct = min(abs(fluct_dn), abs(fluct_up))
                        sym_dn[i] = 1.0 - min_fluct
                        sym_up[i] = 1.0 + min_fluct
                    return sym_dn, sym_up
                
                symmetrized_dn_1, symmetrized_up_1 = symmetrize_binwise(qcd_shape_dn_1, qcd_shape_up_1)
                symmetrized_dn_2, symmetrized_up_2 = symmetrize_binwise(qcd_shape_dn_2, qcd_shape_up_2)
                
                qcd_nuisance_1 = {
                    "pass": self.qcd_pass_1,
                    "loosepass": self.qcd_loosepass_1,
                    "fail": self.qcd_fail_1,
                }[region]

                qcd_nuisance_2 = {
                    "pass": self.qcd_pass_2,
                    "loosepass": self.qcd_loosepass_2,
                    "fail": self.qcd_fail_2,
                }[region]

                # syst_template.setParamEffect(
                #     qcd_nuisance_1,
                #     np.minimum(qcd_shape_dn_1, qcd_shape_up_1),
                #     np.maximum(qcd_shape_dn_1, qcd_shape_up_1),
                # )

                # syst_template.setParamEffect(
                #     qcd_nuisance_2,
                #     np.minimum(qcd_shape_dn_2, qcd_shape_up_2),
                #     np.maximum(qcd_shape_dn_2, qcd_shape_up_2),
                # )
                
                syst_template.setParamEffect(
                    qcd_nuisance_1,
                    qcd_shape_dn_1,
                    qcd_shape_up_1, scale=0.25
                )

                syst_template.setParamEffect(
                    qcd_nuisance_2,
                    qcd_shape_dn_2,
                    qcd_shape_up_2, scale=0.25
                )
                if analysis_region == 'topCR':
                    if self.lowmass_combined:
                        syst_template.setParamEffect(self.qcd_lowmass_top_1[0],
                                                        np.array([symmetrized_dn_1[ix] if ix <= 3 else 1.0 for ix in range(len(qcd['nom']))]),
                                                        np.array([symmetrized_up_1[ix] if ix <= 3 else 1.0 for ix in range(len(qcd['nom']))]))
                        
                        syst_template.setParamEffect(self.qcd_lowmass_top_2[0],
                                                        np.array([symmetrized_dn_2[ix] if ix <= 3 else 1.0 for ix in range(len(qcd['nom']))]),
                                                        np.array([symmetrized_up_2[ix] if ix <= 3 else 1.0 for ix in range(len(qcd['nom']))]))
                    else: 
                        for imx in range(len(self.qcd_lowmass_top_1)):
                            syst_template.setParamEffect(self.qcd_lowmass_top_1[imx],
                                                        np.array([symmetrized_dn_1[imx] if ix == imx else 1.0 for ix in range(len(qcd['nom']))]),
                                                        np.array([symmetrized_up_1[imx] if ix == imx else 1.0 for ix in range(len(qcd['nom']))]))
                        
                        for imx in range(len(self.qcd_lowmass_top_2)):
                            syst_template.setParamEffect(self.qcd_lowmass_top_2[imx],
                                                        np.array([symmetrized_dn_2[imx] if ix == imx else 1.0 for ix in range(len(qcd['nom']))]),
                                                        np.array([symmetrized_up_2[imx] if ix == imx else 1.0 for ix in range(len(qcd['nom']))]))
                elif analysis_region == 'wlnuCR':
                    if self.lowmass_combined:
                        syst_template.setParamEffect(self.qcd_lowmass_wlnu_1[0], 
                                                        np.array([symmetrized_dn_1[ix] if ix <= 3 else 1.0 for ix in range(len(qcd['nom']))]),
                                                        np.array([symmetrized_up_1[ix] if ix <= 3 else 1.0 for ix in range(len(qcd['nom']))]))
                        
                        syst_template.setParamEffect(self.qcd_lowmass_wlnu_2[0],
                                                        np.array([symmetrized_dn_2[ix] if ix <= 3 else 1.0 for ix in range(len(qcd['nom']))]),
                                                        np.array([symmetrized_up_2[ix] if ix <= 3 else 1.0 for ix in range(len(qcd['nom']))]))
                    else:
                        for imx in range(len(self.qcd_lowmass_wlnu_1)):
                            syst_template.setParamEffect(self.qcd_lowmass_wlnu_1[imx], 
                                                        np.array([symmetrized_dn_1[imx] if ix == imx else 1.0 for ix in range(len(qcd['nom']))]),
                                                        np.array([symmetrized_up_1[imx] if ix == imx else 1.0 for ix in range(len(qcd['nom']))]))
                        
                        for imx in range(len(self.qcd_lowmass_wlnu_2)):
                            syst_template.setParamEffect(self.qcd_lowmass_wlnu_2[imx],
                                                        np.array([symmetrized_dn_2[imx] if ix == imx else 1.0 for ix in range(len(qcd['nom']))]),
                                                        np.array([symmetrized_up_2[imx] if ix == imx else 1.0 for ix in range(len(qcd['nom']))]))        
                elif analysis_region == 'SR':
                    if self.lowmass_combined:
                        syst_template.setParamEffect(self.qcd_lowmass_1[0],
                                                        np.array([symmetrized_dn_1[ix] if ix <= 3 else 1.0 for ix in range(len(qcd['nom']))]),
                                                        np.array([symmetrized_up_1[ix] if ix <= 3 else 1.0 for ix in range(len(qcd['nom']))]))

                        syst_template.setParamEffect(self.qcd_lowmass_2[0],
                                                        np.array([symmetrized_dn_2[ix] if ix <= 3 else 1.0 for ix in range(len(qcd['nom']))]),
                                                        np.array([symmetrized_up_2[ix] if ix <= 3 else 1.0 for ix in range(len(qcd['nom']))]))
                    else:
                        for imx in range(len(self.qcd_lowmass_1)):
                            syst_template.setParamEffect(self.qcd_lowmass_1[imx],
                                                    np.array([symmetrized_dn_1[imx] if ix == imx else 1.0 for ix in range(len(qcd['nom']))]),
                                                    np.array([symmetrized_up_1[imx] if ix == imx else 1.0 for ix in range(len(qcd['nom']))]))

                        for imx in range(len(self.qcd_lowmass_2)):
                            syst_template.setParamEffect(self.qcd_lowmass_2[imx],
                                                    np.array([symmetrized_dn_2[imx] if ix == imx else 1.0 for ix in range(len(qcd['nom']))]),
                                                    np.array([symmetrized_up_2[imx] if ix == imx else 1.0 for ix in range(len(qcd['nom']))])) 


                # for imx in range(len(self.qcd_lowmass_top)):
                #     syst_template.setParamEffect(self.qcd_lowmass_top[imx],
                #                         np.array([qcd_shape_dn_1[imx] if ix==imx else 1.0 for ix in range(len(qcd['nom']))]),
                #                         np.array([qcd_shape_up_1[imx] if ix==imx else 1.0 for ix in range(len(qcd['nom']))]))
            

                # for imx in range(len(self.qcd_lowmass_wlnu)):
                #     syst_template.setParamEffect(self.qcd_lowmass_wlnu[imx],
                #                             np.array([qcd_shape_dn_1[imx] if ix==imx else 1.0 for ix in range(len(qcd['nom']))]),
                #                             np.array([qcd_shape_up_1[imx] if ix==imx else 1.0 for ix in range(len(qcd['nom']))]))
                
                '''
                for imx in range(len(self.qcd_lowmass_1)):
                    syst_template.setParamEffect(self.qcd_lowmass_1[imx],
                                            np.array([qcd_shape_dn_1[imx] if ix==imx else 1.0 for ix in range(len(qcd['nom']))]),
                                            np.array([qcd_shape_up_1[imx] if ix==imx else 1.0 for ix in range(len(qcd['nom']))]))
                
                for imx in range(len(self.qcd_lowmass_2)):
                    syst_template.setParamEffect(self.qcd_lowmass_2[imx],
                                            np.array([qcd_shape_dn_2[imx] if ix==imx else 1.0 for ix in range(len(qcd['nom']))]),
                                            np.array([qcd_shape_up_2[imx] if ix==imx else 1.0 for ix in range(len(qcd['nom']))]))
                '''
                # bin by bin nuisances
                #  for imx in range(len(self.qcd_lowmass_1)):
                #     syst_template.setParamEffect(self.qcd_lowmass_1[imx],
                #                                 np.array([symmetrized_dn_1[imx] if ix == imx else 1.0 for ix in range(len(qcd['nom']))]),
                #                                 np.array([symmetrized_up_1[imx] if ix == imx else 1.0 for ix in range(len(qcd['nom']))]))

                # for imx in range(len(self.qcd_lowmass_2)):
                #     syst_template.setParamEffect(self.qcd_lowmass_2[imx],
                #                                 np.array([symmetrized_dn_2[imx] if ix == imx else 1.0 for ix in range(len(qcd['nom']))]),
                #                                 np.array([symmetrized_up_2[imx] if ix == imx else 1.0 for ix in range(len(qcd['nom']))]))

                # syst_template.setParamEffect(self.qcd_lowmass_top, qcd_shape_dn, qcd_shape_up)
                    
                # syst_template.setParamEffect(self.qcd_lowmass_wlnu, qcd_shape_dn, qcd_shape_up)

        return syst_template

    def build_channel(
        self,
        analysis_region,
        region,
        hchannel,
        qcd,
        unblind=False,
        singlebin=False,
        rebin=False,
        qcdnormSF=None,
        debug=False,
    ):
        cat = self.__name
        
        analysis_region_str = "" if analysis_region == "SR" else analysis_region

        ch = rl.Channel(f"{analysis_region_str}{region}{cat}{self.year}")

        vals = {}

        def smass(sName):
            if 'h125' in sName or 'htt125' in sName:
                _mass = 125.
            elif sName in ['wjets', 'vv', 'vqq', 'wlnu']:
                _mass = 80.379
            elif sName in ['ztt', 'zem', 'dy']:
                _mass = 91.
            else:
                raise ValueError("DAFUQ is {}".format(sName))
            return _mass

        for sample in hchannel.identifiers("sample"):
            if sample.name == "ignore":
                continue
            if (
                sample.name == "data_obs"
                and region == "pass"
                and not unblind
                and analysis_region == "SR"
            ):
                continue

            logging.debug(
                f"Building template for sample {sample} and region {region} and analysis_region {analysis_region}, singlebin {singlebin}"
            )
            print(f"Building template for sample {sample} and region {region} and analysis_region {analysis_region}, singlebin {singlebin}")
            # get template from that MC/data sample
            h = hchannel[sample]
            template, events = self.get_template(
                h, region, "nominal", singlebin, rebin, debug, clip=True
            )

            # if the sample is multijet(QDC) take "qcd" i.e. the prediction from data
            if sample.name == "multijet":
                clip = (
                    np.clip(np.array([np.sum(qcd["nom"])]), 0.0, None)
                    if singlebin
                    else np.clip(qcd["nom"], 0.0, None)
                )
    
                template = (
                    clip,
                    (self.mttone.binning if singlebin else (self.mttrebin.binning if rebin else self.mtt.binning)),
                    (self.mttone.name if singlebin else (self.mttrebin.name if rebin else self.mtt.name)),
                    (np.array([template[-1]]) if singlebin else template[-1])
                )

            if sample.name == "data_obs":
                #print("XXXX")
                #print("Data", template)
                print(template)
                ch.setObservation(template[:-1])
                #ch.setObservation(template, read_sumw2=True)
            else:
                print(template)
                print(sample.name)
                sample_template = rl.TemplateSample(
                    f"{ch.name}_{sample.name}",
                    rl.Sample.SIGNAL
                    if self.signame in sample.name
                    else rl.Sample.BACKGROUND,
                    template,
                )

            
                # MAYBE add in automcstats soon? 
                # if 'multijet' not in sample.name:
                #     sample_template.autoMCStats(epsilon=1e-4, lnN=False)


                # shape systematics
                sample_template = self.systs_shape(
                    h, region, sample, sample_template, singlebin
                )
                # mass systematics - skipping for now
                sample_template = self.systs_mass(
                    h, region, sample, sample_template, template, singlebin, analysis_region, 
                )
                # norm systematics
                sample_template = self.systs_norm(
                    region, events, qcd, sample, sample_template, singlebin, analysis_region
                )
                # QCD norm SF
                if sample.name == "multijet" and qcdnormSF is not None:
                    sample_template.setParamEffect(qcdnormSF, 1 * qcdnormSF)

                # add sample
                ch.addSample(sample_template)

        # blind signal region
        if not unblind and region == "pass" and analysis_region == "SR":
            ch.setObservation(
                (
                    np.zeros(len(self.mttone.binning) - 1)
                    if singlebin
                    else np.zeros(len(self.mtt.binning) - 1),
                    self.mttone.binning if singlebin else self.mtt.binning,
                    self.mttone.name if singlebin else self.mtt.name,
                )
            )
        ch.autoMCStats(epsilon=1e-4, threshold=1)
        #ch.autoMCStats(epsilon=1e-4)


        #masking 
        # if analysis_region == 'wlnuCR' and region == "fail":
        #     mask = np.array([False, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True])
        #     ch.mask = mask
        
        # if analysis_region == 'wlnuCR' and region == "pass":
        #     mask = np.array([False, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True])
        #     ch.mask = mask

        # if analysis_region == 'topCR' and region == "pass":
        #     mask = np.array([True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True])
        #     ch.mask = mask
    
        # add channel
        self.model.addChannel(ch)
        
        return vals

    def set_expectation(self, str_fail, str_loose, str_pass):
        
        # TODO: Introduce unifiedBkgEff if necessary?

        # top normalization
        # QUESTION: shouldn't this be loose/pass?
        unifiedBkgEff_flag = False

        topLPF = (
            self.model[str_pass]["top"].getExpectation(nominal=True).sum()
            / self.model[str_loose]["top"].getExpectation(nominal=True).sum()
        )
        topPF = (
            self.model[str_loose]["top"].getExpectation(nominal=True).sum()
            + self.model[str_pass]["top"].getExpectation(nominal=True).sum()
        ) / self.model[str_fail]["top"].getExpectation(nominal=True).sum()

        if unifiedBkgEff_flag: 

            self.model[str_pass]["top"].setParamEffect(self.bkgLeffSF, 1 * self.bkgLeffSF)
            self.model[str_loose]["top"].setParamEffect(
                self.bkgLeffSF, (1 - self.bkgLeffSF) * topLPF + 1
            )


            self.model[str_loose]["top"].setParamEffect(self.bkgeffSF, 1 * self.bkgeffSF)
            self.model[str_pass]["top"].setParamEffect(self.bkgeffSF, 1 * self.bkgeffSF)
            self.model[str_fail]["top"].setParamEffect(
                self.bkgeffSF, (1 - self.bkgeffSF) * topPF + 1
            )

        else: 
            self.model[str_pass]["top"].setParamEffect(self.topLeffSF, 1 * self.topLeffSF)
            self.model[str_loose]["top"].setParamEffect(
                self.topLeffSF, (1 - self.topLeffSF) * topLPF + 1
            )

            self.model[str_loose]["top"].setParamEffect(self.topeffSF, 1 * self.topeffSF)
            self.model[str_pass]["top"].setParamEffect(self.topeffSF, 1 * self.topeffSF)
            self.model[str_fail]["top"].setParamEffect(
                self.topeffSF, (1 - self.topeffSF) * topPF + 1
            )

        # wlnu normalization
        wlnuPF = (
            self.model[str_loose]["wlnu"].getExpectation(nominal=True).sum()
            + self.model[str_pass]["wlnu"].getExpectation(nominal=True).sum()
        ) / self.model[str_fail]["wlnu"].getExpectation(nominal=True).sum()
        wlnuLPF = (
            self.model[str_pass]["wlnu"].getExpectation(nominal=True).sum()
            / self.model[str_loose]["wlnu"].getExpectation(nominal=True).sum()
        )
        # pass/(loose + pass)
        wlnuRLPF = 1.0 / (1.0 + (1.0 / wlnuLPF))

        if unifiedBkgEff_flag: 
            self.model[str_pass]["wlnu"].setParamEffect(
            self.bkgLeffSF, 1 * self.bkgLeffSF
            )
            self.model[str_loose]["wlnu"].setParamEffect(
                self.bkgLeffSF, (1 - self.bkgLeffSF) * wlnuLPF + 1
            )
            self.model[str_loose]["wlnu"].setParamEffect(self.bkgeffSF, 1 * self.bkgeffSF)
            self.model[str_pass]["wlnu"].setParamEffect(self.bkgeffSF, 1 * self.bkgeffSF)
            self.model[str_fail]["wlnu"].setParamEffect(
                self.bkgeffSF, (1 - self.bkgeffSF) * wlnuPF + 1
            )

        else:

            self.model[str_pass]["wlnu"].setParamEffect(
            self.wlnuLeffSF, 1 * self.wlnuLeffSF
            )
            self.model[str_loose]["wlnu"].setParamEffect(
                self.wlnuLeffSF, (1 - self.wlnuLeffSF) * wlnuLPF + 1
            )

            self.model[str_loose]["wlnu"].setParamEffect(self.wlnueffSF, 1 * self.wlnueffSF)
            self.model[str_pass]["wlnu"].setParamEffect(self.wlnueffSF, 1 * self.wlnueffSF)
            self.model[str_fail]["wlnu"].setParamEffect(
                self.wlnueffSF, (1 - self.wlnueffSF) * wlnuPF + 1
            )

        # dy normalization
        dyLP = (
            self.model[str_pass]["dy"].getExpectation(nominal=True).sum()
            / self.model[str_loose]["dy"].getExpectation(nominal=True).sum()
        )
        self.model[str_pass]["dy"].setParamEffect(self.dy_eff, 1 * self.dy_eff)
        self.model[str_loose]["dy"].setParamEffect(
            self.dy_eff, (1 - self.dy_eff) * dyLP + 1
        )

        # htt
        httLP = (
            self.model[str_pass]["htt125"].getExpectation(nominal=True).sum()
            / self.model[str_loose]["htt125"].getExpectation(nominal=True).sum()
        )
        self.model[str_pass]["htt125"].setParamEffect(self.dy_eff, 1 * self.dy_eff)
        # self.model[str_loose]["htt125"].setParamEffect(
        #     self.dy_eff, (1 - self.dy_eff) * httLP + 1
        # )
        # phitt
        if not self.doHtt:
            for m in self.masspoints:
                phittLP = (
                    self.model[str_pass][f"phitt{m}"].getExpectation(nominal=True).sum()
                    / self.model[str_loose][f"phitt{m}"]
                    .getExpectation(nominal=True)
                    .sum()
                    if self.model[str_loose][f"phitt{m}"]
                    .getExpectation(nominal=True)
                    .sum()
                    > 0.0
                    else 1.0
                )
                self.model[str_pass][f"phitt{m}"].setParamEffect(
                    self.dy_eff, 1 * self.dy_eff
                )
                self.model[str_loose][f"phitt{m}"].setParamEffect(
                    self.dy_eff, (1 - self.dy_eff) * phittLP + 1
                )
            self.model[str_pass]["htt125"].setParamEffect(self.rh125, 1.10)
            self.model[str_loose]["htt125"].setParamEffect(self.rh125, 1.10)
            self.model[str_fail]["htt125"].setParamEffect(self.rh125, 1.10)
