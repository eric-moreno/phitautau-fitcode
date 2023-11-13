import numpy as np
import logging
import rhalphalib as rl
from utils import intRegion, getQCDFromData
from coffea import hist

rl.ParametericSample.PreferRooParametricHist = True


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

        # bin slice
        self.lowbin = 2
        self.highbin = -1
        self.mttbins_nom = self.mttbins[self.lowbin : self.highbin]
        self.mttrange = range(len(self.mttbins_nom))

        # mass range
        self.mass = "massreg"
        self.lowqcdmass = 55.0 if self.islephad else 105.0
        self.highmass = 145.0 if self.islephad else 145.0
        self.lowqcdincrease = 0.5
        self.highbkgincrease = 0.3

        # pt range
        self.pt_min = 300

        # signal names
        self.doHtt = True
        self.signame = "htt125" if self.doHtt else "phitt"
        self.sigexname = "htt125" if self.doHtt else "phitt50"

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

        # group samples
        self.sample_groups = {
            "data_obs": ["data"],
            "top": ["tt-had", "tt-semilep", "tt-dilep", "st"],
            "htt125": ["h125"],
            "multijet": ["qcd"],
            "dy": ["zem", "ztt"],
            "wlnu": ["wjets", "vv", "vqq"],
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

        # met cuts
        cuts = {
            "lowmet": 20.0,
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
                    "sig": slices["lowtohighmet"],
                    "top_cr": slices["lowtohighmet"],
                    "wlnu_cr": slices["lowmet"],
                    "qcd_cr": slices["lowmet"],
                },
                "fail": {
                    "sig": slices["lowtohighmet"],
                    "top_cr": slices["lowtohighmet"],
                    "wlnu_cr": slices["lowmet"],
                    "qcd_cr": slices["lowmet"],
                },
            }

        # nn regions
        self.nnregions = ["fail", "loosepass", "pass"]

        # variation of background yields for systematics
        self.variations = ["nom", "dn", "up"]

        # initialize model
        self.model = rl.Model(f"{cat}Model")

        # initialize observable
        self.mtt = rl.Observable(self.mass, self.mttbins_nom)
        self.mttone = rl.Observable(
            self.mass,
            np.array([self.mttbins[self.lowbin], self.mttbins[self.highbin - 1]]),
        )

        # initialize nuisances
        self.nuisances(no_syst)

    def nuisances(self, no_syst: bool):
        self.logger.info(f"Initializing nuisances")
        cat = self.__name

        # add UE and JES/JER
        syst_dict = {
            samp: {
                "uescale": [
                    rl.NuisanceParameter("uescale", "shape"),
                    "UESDown",
                    "UESUp",
                ],
                "jescale": [
                    rl.NuisanceParameter("jescale", "shape"),
                    "JESDown",
                    "JESUp",
                ],
                "jeresol": [
                    rl.NuisanceParameter("jeresol", "shape"),
                    "JERDown",
                    "JERUp",
                ],
            }
            for samp in self.sample_groups
            if samp not in ["data_obs", "ignore", "multijet"]
        }

        # add toppt uncertainty for top region
        syst_dict["top"]["toppt"] = [
            rl.NuisanceParameter("toppt", "shape"),
            "nominal",
            "TopPtReweightUp",
        ]

        # add prefire uncertainty
        if self.year not in ["2018"]:
            for sample in syst_dict:
                syst_dict[sample]["l1prefire"] = [
                    rl.NuisanceParameter("l1prefire", "shape"),
                    "L1PreFiringDown",
                    "L1PreFiringUp",
                ]

        # dictionary of systematic uncertainty variations to build from existing mass templates
        if no_syst:
            self.syst_dict = {}
        else:
            self.syst_dict = syst_dict

        # QCD normalization
        self.qcdnormSF = rl.IndependentParameter(f"qcdnormSF_{cat}", 1.0, 0, 10)
        self.qcdnormSF_top = rl.IndependentParameter(f"qcdnormSF_top_{cat}", 1.0, 0, 10)
        self.qcdnormSF_wlnu = rl.IndependentParameter(
            f"qcdnormSF_wlnu_{cat}", 1.0, 0, 10
        )

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
            rl.NuisanceParameter(f"wlnuhighmass_bin{ix}_{cat}", "shape")
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
        self.qcd_fail = rl.NuisanceParameter(f"qcd_Rfail_{cat}", "shape")
        self.qcd_loosepass = rl.NuisanceParameter(f"qcd_Rloosepass_{cat}", "shape")
        self.qcd_pass = rl.NuisanceParameter(f"qcd_Rpass_{cat}", "shape")

        # add uncertainty for low mass
        self.qcd_lowmass = [
            rl.NuisanceParameter(f"qcd_lowmass_bin{ix}_{cat}", "shape")
            for ix in self.mttrange
            if self.mttbins_nom[ix] < self.lowqcdmass
        ]
        self.qcd_lowmass_top = [
            rl.NuisanceParameter(f"qcd_lowmass_top_bin{ix}_{cat}", "shape")
            for ix in self.mttrange
            if self.mttbins_nom[ix] < self.lowqcdmass
        ]
        self.qcd_lowmass_wlnu = [
            rl.NuisanceParameter(f"qcd_lowmass_wlnu_bin{ix}_{cat}", "shape")
            for ix in self.mttrange
            if self.mttbins_nom[ix] < self.lowqcdmass
        ]

    def get_qcd(self, hists: dict, region: str, default: float = 1.0):
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
            hists, region, self.nnCut, self.nnCutLoose, default=default
        )

    def _get_region(
        self, h: hist.Hist, region: str, syst: str = "nominal", debug: bool = False
    ):
        """
        Get region by integrating over nn cut and/or met/mass/pt cuts
        """
        return intRegion(
            h,
            region,
            self.nnCut,
            self.nnCutLoose,
            systematic=syst,
            debug=debug,
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

    def get_template(self, h, region, syst, singlebin, debug, clip=True):
        """
        Get template ntuple
        """
        tempint = self._get_region(h, region, syst, debug)
        events = self._events(tempint, singlebin, clip)
        sumw = self._sumw(tempint, singlebin)
        template = (
            events,
            self.mttone.binning if singlebin else self.mtt.binning,
            self.mttone.name if singlebin else self.mtt.name,
            sumw,
        )
        return template, events

    def systs_mass(self, h, region, sample, sample_template, singlebin=False):
        """
        Add mass systematics
        """
        cat = self.__name
        syst_template = sample_template
        if (
            sample.name == "htt125" or sample.name == "dy" or "phi" in sample.name
        ) and not singlebin:
            # why not use rl mass shift functions?
            m_scale = rl.NuisanceParameter(f"massscale_{cat}", "shape")

            nominal = self._get_region(h, region)
            nom = self._events(nominal, clip=False)
            nom_full = nominal[0]
            shift_dn = nom_full[self.lowbin - 1 : self.highbin - 1]
            shift_up = nom_full[
                self.lowbin + 1 : self.highbin + 1 if self.highbin != -1 else None
            ]
            shiftwidth = self.mttbins[1] - self.mttbins[0]

            dnfrac = np.array(
                [
                    shiftwidth
                    / (
                        self.mttbins[self.lowbin + ib]
                        - self.mttbins[self.lowbin + ib - 1]
                    )
                    for ib in range(len(nom))
                ]
            )
            shift_dn = shift_dn * dnfrac + nom * (1.0 - dnfrac)

            upfrac = np.array(
                [
                    shiftwidth
                    / (
                        self.mttbins[self.lowbin + ib + 2]
                        - self.mttbins[self.lowbin + ib + 1]
                    )
                    for ib in range(len(nom))
                ]
            )
            shift_up = shift_up * upfrac + nom * (1.0 - upfrac)

            syst_template.setParamEffect(
                m_scale,
                np.divide(shift_dn, nom, out=np.ones_like(nom), where=nom > 0.0),
                np.divide(shift_up, nom, out=np.ones_like(nom), where=nom > 0.0),
            )

        return syst_template

    def systs_shape(self, h, region, sample, sample_template, singlebin=False):
        """
        Add systematics with shape samples
        """

        syst_template = sample_template

        if sample.name in self.syst_dict:
            tempint = self._get_region(h, region)
            nom = self._events(tempint)

            for syst in self.syst_dict[sample.name]:
                nuisance, syst_dn, syst_up = self.syst_dict[sample.name][syst]
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
        return syst_template

    def systs_norm(self, region, events, qcd, sample, sample_template, singlebin=False):
        """
        Add systematics for normalization
        """

        cat = self.__name
        syst_template = sample_template

        if sample.name == "top":
            syst_template.setParamEffect(
                rl.NuisanceParameter("CMS_top_norm", "lnN"), 1.05
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
                rl.NuisanceParameter("CMS_wlnu_norm", "lnN"), 1.10
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

        if sample.name == "vvqq":
            syst_template.setParamEffect(
                rl.NuisanceParameter("CMS_vvqq_norm", "lnN"), 1.20
            )
        if sample.name == "dy":
            syst_template.setParamEffect(
                rl.NuisanceParameter("CMS_vvqq_norm", "lnN"), 1.05
            )

        lumi_16 = rl.NuisanceParameter("CMS_lumi_16", "lnN")
        lumi_17 = rl.NuisanceParameter("CMS_lumi_17", "lnN")
        lumi_18 = rl.NuisanceParameter("CMS_lumi_18", "lnN")
        lumi_all = rl.NuisanceParameter("CMS_lumi_all", "lnN")
        lumi_1718 = rl.NuisanceParameter("CMS_lumi_1718", "lnN")
        lumi_list = {
            "2016": [lumi_16, lumi_all],
            "2017": [lumi_17, lumi_all, lumi_1718],
            "2018": [lumi_18, lumi_all, lumi_1718],
        }[self.year]
        lumi_vals = {
            "2016": [1.01, 1.006],
            "2017": [1.02, 1.009, 1.006],
            "2018": [1.015, 1.02, 1.002],
        }[self.year]

        if sample.name != "multijet":
            syst_template.setParamEffect(
                rl.NuisanceParameter(f"CMS_trig_{cat}", "lnN"), 1.02
            )
            if self.islephad:
                syst_template.setParamEffect(
                    rl.NuisanceParameter(f"CMS_id_{cat}", "lnN"), 1.02
                )
            for il, lumi in enumerate(lumi_list):
                syst_template.setParamEffect(lumi, lumi_vals[il])
        else:
            if not singlebin:
                qcd_shape_dn = np.divide(
                    qcd["dn"],
                    qcd["nom"],
                    out=np.ones_like(qcd["nom"]),
                    where=qcd["nom"] > 0.0,
                )
                qcd_shape_up = np.divide(
                    qcd["up"],
                    qcd["nom"],
                    out=np.ones_like(qcd["nom"]),
                    where=qcd["nom"] > 0.0,
                )
                qcd_nuisance = {
                    "pass": self.qcd_pass,
                    "loosepass": self.qcd_loosepass,
                    "fail": self.qcd_fail,
                }[region]
                syst_template.setParamEffect(
                    qcd_nuisance,
                    np.minimum(qcd_shape_dn, qcd_shape_up),
                    np.maximum(qcd_shape_dn, qcd_shape_up),
                )
                # for imx in range(len(qcd_lowmass_wlnu)):
                #    sample.setParamEffect(qcd_lowmass_wlnu[imx],
                #                          np.array([qcd_shape_dn[imx] if ix==imx else 1. for ix in range(len(qcdpred))]),
                #                          np.array([qcd_shape_up[imx] if ix==imx else 1. for ix in range(len(qcdpred))]))

        return syst_template

    def build_channel(
        self,
        analysis_region,
        region,
        hchannel,
        qcd,
        unblind=False,
        singlebin=False,
        qcdnormSF=None,
        debug=False,
    ):
        cat = self.__name
        
        analysis_region_str = "" if analysis_region == "SR" else analysis_region

        ch = rl.Channel(f"{analysis_region_str}{region}{cat}{self.year}")

        vals = {}
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
            # get template from that MC/data sample
            h = hchannel[sample]
            template, events = self.get_template(
                h, region, "nominal", singlebin, debug, clip=True
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
                    self.mttone.binning if singlebin else self.mtt.binning,
                    self.mttone.name if singlebin else self.mtt.name,
                    np.array([1000000.0])
                    if singlebin
                    else np.ones_like(qcd["nom"]) * 1000000.0,
                )

            if sample.name == "data_obs":
                ch.setObservation(template, read_sumw2=True)
            else:
                sample_template = rl.TemplateSample(
                    f"{ch.name}_{sample.name}",
                    rl.Sample.SIGNAL
                    if self.signame in sample.name
                    else rl.Sample.BACKGROUND,
                    template,
                )
                # shape systematics
                sample_template = self.systs_shape(
                    h, region, sample, sample_template, singlebin
                )
                # mass systematics - skipping for now
                # sample_template = self.systs_mass(
                #    h, region, ptslice, sample, sample_template, singlebin
                # )
                # norm systematics
                sample_template = self.systs_norm(
                    region, events, qcd, sample, sample_template, singlebin
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

        # add channel
        self.model.addChannel(ch)
        
        return vals

    def set_expectation(self, str_fail, str_loose, str_pass):
        print(str_pass)
        # TODO: Introduce unifiedBkgEff if necessary?

        # top normalization
        # QUESTION: shouldn't this be loose/pass?
        topLPF = (
            self.model[str_pass]["top"].getExpectation(nominal=True).sum()
            / self.model[str_loose]["top"].getExpectation(nominal=True).sum()
        )
        topPF = (
            self.model[str_loose]["top"].getExpectation(nominal=True).sum()
            + self.model[str_pass]["top"].getExpectation(nominal=True).sum()
        ) / self.model[str_fail]["top"].getExpectation(nominal=True).sum()

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
        self.model[str_loose]["htt125"].setParamEffect(
            self.dy_eff, (1 - self.dy_eff) * httLP + 1
        )
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
            self.model[str_loose]["htt125"].setParamEffect(self.rh125, 1.10)
