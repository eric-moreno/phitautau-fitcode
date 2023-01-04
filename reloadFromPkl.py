from __future__ import print_function, division
import sys
import os
import rhalphalib as rl
import numpy as np
import scipy.stats
import argparse
import pickle
import ROOT
rl.util.install_roofit_helpers()
#rl.ParametericSample.PreferRooParametricHist = False
rl.ParametericSample.PreferRooParametricHist = True

def reloadPkl(args):

    outdirname = args.outdir[-1]

    if len(args.outdir)>1:
        args.outdir = args.outdir[:-1]
        singleDir = False
    else:
        args.outdir = args.outdir[0]
        singleDir = True

    os.system('mkdir -p %s'%outdirname)
    if singleDir:
        os.system("cp %s/%s/info.txt %s/"%(args.indir,args.outdir,args.outdir))
    for m in args.pkl:
        chanlist = []

        if m!='fullModel':
            if singleDir:
                os.system("cp %s/%s/%s.pkl %s/"%(args.indir,args.outdir,m,outdirname))
            else:
                for od in args.outdir:
                    os.system("cp %s/%s/%s.pkl %s/%s%s.pkl"%(args.indir,od,m,outdirname,od.replace('/',''),m))

        for mass in args.masses:
            if m!='fullModel':
                if singleDir:
                    model = pickle.load(open("%s/%s/%s.pkl"%(args.indir,args.outdir,m), "rb"))
                else:
                    model = pickle.load(open("%s/%s/%s.pkl"%(args.indir,args.outdir[0],m), "rb"))
                    for iod in range(1,len(args.outdir)):
                        othermodel = pickle.load(open("%s/%s/%s.pkl"%(args.indir,args.outdir[iod],m), "rb"))
                        for c in othermodel:
                            model.addChannel(othermodel[c.name])
            else:
                if singleDir:
                    model = pickle.load(open("%s/%s/%s.pkl"%(args.indir,args.outdir,'hadelModel'), "rb"))
                    for om in ['hadmuModel','hadhadModel']:
                        othermodel = pickle.load(open("%s/%s/%s.pkl"%(args.indir,args.outdir,om), "rb"))
                        for c in othermodel:
                            model.addChannel(othermodel[c.name])
                else:
                    model = pickle.load(open("%s/%s/%s.pkl"%(args.indir,args.outdir[0],'hadelModel'), "rb"))
                    isFirst = True
                    for iod in range(0,len(args.outdir)):
                        for om in ['hadelModel','hadmuModel','hadhadModel']:
                            if isFirst:
                                isFirst = False
                                continue
                            othermodel = pickle.load(open("%s/%s/%s.pkl"%(args.indir,args.outdir[iod],om), "rb"))
                            for c in othermodel:
                                model.addChannel(othermodel[c.name])
                model._name = 'fullModel'
            for c in model:
                chanlist.append(c.name)
                remlist = []
                for s in model[c.name]:
                    if 'phitt' in s.name and not s.name.endswith('phitt%s'%mass):
                        remlist.append(s)
                for s in remlist:
                    del model[c.name]._samples[s.name]
                #model[c.name].autoMCStats(channel_name=c.name)

            model.renderCombine("%s/%s_m%s/"%(outdirname,m,mass))
    
            if args.verbose:
                year = args.year
                print('top',model['ptbin0pass%s'%year]['top'].getExpectation(nominal=True))
                print('wlnu',model['ptbin0pass%s'%year]['wlnu'].getExpectation(nominal=True))
                print('ztt',model['ptbin0pass%s'%year]['ztt'].getExpectation(nominal=True))
                print('htt125',model['ptbin0pass%s'%year]['htt125'].getExpectation(nominal=True))
                if not args.dohtt:
                    print('phitt%s'%mass,model['ptbin0pass%s'%year]['phitt%s'%mass].getExpectation(nominal=True))
    
            #mypath = "%s/%s"%(args.outdir,m)
            #fullfiles = ["%s/%s"%(mypath,f) for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath, f)) and 'txt' in f]
            #for f in fullfiles:
            #    os.system("echo \* autoMCStats 0 >> %s"%f)

            #os.system('sed -i "1a\\echo \'* autoMCStats 0\' >> %s_combined.txt" %s/%s_m%s/build.sh'%(m,args.outdir,m,mass))
            #os.system("echo combine -M %s %s_combined.root -m %s --rMin -1 --rMax 5 -t -1 --expectSignal %s -n .%s -v 3 %s >> %s/%s/build.sh"%(args.mode,m,mass,args.expsig,m,args.addargs,args.outdir,m))
            if '2016APV' in args.outdir:
                year = '2016APV'
            elif '2016' in args.outdir:
                year = '2016'
            elif '2017' in args.outdir:
                year = '2017'
            elif '2018' in args.outdir:
                year = '2018'

            os.system("echo '' >> %s/%s_m%s/build.sh"%(outdirname,m,mass))
            os.system("cp %s/%s_m%s/build{,_gof}.sh"%(outdirname,m,mass))
            os.system("sed -i 's/text2workspace.py.*$/& --channel-masks/' %s/%s_m%s/build_gof.sh"%(args.outdir,m,mass))

            if args.nopass:
                os.system("sed -i 's/ptbin[0-9]passhad[muelhad]*%s=ptbin[0-9]passhad[muelhad]*%s.txt//g'  %s/%s_m%s/build.sh"%(year,year,outdirname,m,mass))

            maskstr = ''
            if args.mask:
                os.system("sed -i 's/text2workspace.py.*$/& --channel-masks/' %s/%s_m%s/build.sh"%(args.outdir,m,mass))
                maskstr = '--setParameters mask_ptbin0passhadel{0}=1,mask_ptbin1passhadel{0}=1,mask_ptbin0passhadmu{0}=1,mask_ptbin1passhadmu{0}=1,mask_ptbin0passhadhad{0}=1,mask_ptbin1passhadhad{0}=1,mask_ptbin2passhadhad{0}=1,mask_ptbin3passhadhad{0}=1'.format(year)
            #os.system("echo combine -M %s %s_combined.root -m %s --rMin -1 --rMax 5 -t -1 -n .%s -v 3 %s %s >> %s/%s_m%s/build.sh"%(args.mode,m,mass,m,args.addargs,maskstr,args.outdir,m,mass))
            #if args.full:
            #    os.system("echo combine -M FitDiagnostics %s_combined.root -m %s --rMin -1 --rMax 5 -t -1 --expectSignal %s --saveShapes --saveWithUncertainties -n .%s -v 3 %s %s >> %s/%s_m%s/build.sh"%(m,mass,args.expsig,m,args.addargs,maskstr,args.outdir,m,mass))
            #    os.system("echo combineTool.py -M Impacts -d %s_combined.root -m %s --rMin -1 --rMax 5 --robustFit 1 --doInitialFit --expectSignal %s -t -1 %s %s >> %s/%s_m%s/build.sh"%(m,mass,args.expsig,args.addargs,maskstr,args.outdir,m,mass))
            #    os.system("echo combineTool.py -M Impacts -d %s_combined.root -m %s --rMin -1 --rMax 5 --robustFit 1 --doFits --expectSignal %s -t -1 %s %s >> %s/%s_m%s/build.sh"%(m,mass,args.expsig,args.addargs,maskstr,args.outdir,m,mass))
            #    os.system("echo combineTool.py -M Impacts -d %s_combined.root -m %s --rMin -1 --rMax 5 --robustFit 1 --output impacts_m%s.json --expectSignal %s -t -1 %s %s >> %s/%s_m%s/build.sh"%(m,mass,mass,args.expsig,args.addargs,maskstr,args.outdir,m,mass))

            os.system("echo combine -M %s %s_combined.root -m %s --rMin -200 --rMax 1000 --expectSignal %s -n .%s -v 3 %s %s >> %s/%s_m%s/build.sh"%(args.mode,m,mass,args.expsig,m,args.addargs,maskstr,outdirname,m,mass))

            if args.full:
                os.system("echo combine -M FitDiagnostics %s_combined.root -m %s --rMin -200 --rMax 1000 --expectSignal %s --saveShapes --saveWithUncertainties -n .%s -v 3 %s %s >> %s/%s_m%s/build.sh"%(m,mass,args.expsig,m,args.addargs,maskstr,outdirname,m,mass))
                os.system("echo combineTool.py -M Impacts -d %s_combined.root -m %s --rMin -200 --rMax 1000 --robustFit 1 --doInitialFit --expectSignal %s %s %s >> %s/%s_m%s/build.sh"%(m,mass,args.expsig,args.addargs,maskstr,outdirname,m,mass))
                os.system("echo combineTool.py -M Impacts -d %s_combined.root -m %s --rMin -200 --rMax 1000 --robustFit 1 --doFits --expectSignal %s %s %s >> %s/%s_m%s/build.sh"%(m,mass,args.expsig,args.addargs,maskstr,outdirname,m,mass))
                os.system("echo combineTool.py -M Impacts -d %s_combined.root -m %s --rMin -200 --rMax 1000 --robustFit 1 --output impacts_m%s.json --expectSignal %s %s %s >> %s/%s_m%s/build.sh"%(m,mass,mass,args.expsig,args.addargs,maskstr,outdirname,m,mass))
                os.system("echo plotImpacts.py -i impacts_m%s.json -o impacts_%s >> %s/%s_m%s/build.sh"%(mass,m,outdirname,m,mass))
                #need to use \# not # for making comments

            os.system("echo combine -M GoodnessOfFit -d %s_combined.root -m %s --algo=saturated -n _result_bonly_CRonly --setParametersForFit mask_ptbin0passhadel%s=1 --setParametersForEval mask_ptbin0passhadel%s=0 --freezeParameters r --setParameters r=0,mask_ptbin0passhadel%s=1 -t 500 --toysFrequentist >> %s/%s_m%s/build_gof.sh"%(m,mass,year,year,year,outdirname,m,mass))
            os.system("echo combine -M GoodnessOfFit -d %s_combined.root -m %s --algo=saturated -n _result_bonly_CRonly_data --setParametersForFit mask_ptbin0passhadel%s=1 --setParametersForEval mask_ptbin0passhadel%s=0 --freezeParameters r --setParameters r=0,mask_ptbin0passhadel%s=1 >> %s/%s_m%s/build_gof.sh"%(m,mass,year,year,year,outdirname,m,mass))
            os.system("echo combine -M GoodnessOfFit -d %s_combined.root -m %s --algo=saturated -n _result_sb -t 500 --toysFrequentist >> %s/%s_m%s/build_gof.sh"%(m,mass,outdirname,m,mass))
            os.system("echo combine -M GoodnessOfFit -d %s_combined.root -m %s --algo=saturated -n _result_sb_data >> %s/%s_m%s/build_gof.sh"%(m,mass,outdirname,m,mass))
            os.system("echo combineTool.py -M CollectGoodnessOfFit --input higgsCombine_result_bonly_CRonly_data.GoodnessOfFit.mH%s.root higgsCombine_result_bonly_CRonly.GoodnessOfFit.mH%s.123456.root -m %s -o gof_CRonly.json >> %s/%s_m%s/build_gof.sh"%(mass,mass,mass,outdirname,m,mass))
            os.system("echo combineTool.py -M CollectGoodnessOfFit --input higgsCombine_result_sb_data.GoodnessOfFit.mH%s.root higgsCombine_result_sb.GoodnessOfFit.mH%s.123456.root -m %s -o gof.json >> %s/%s_m%s/build_gof.sh"%(mass,mass,mass,outdirname,m,mass))
            os.system("echo plotGof.py gof_CRonly.json --statistic saturated --mass %s.0 -o gof_plot_CRonly --title-right='%s_%s_CRonly' >> %s/%s_m%s/build_gof.sh"%(mass,m,year,outdirname,m,mass))
            os.system("echo plotGof.py gof.json --statistic saturated --mass %s.0 -o gof_plot --title-right='%s_%s' >> %s/%s_m%s/build_gof.sh"%(mass,m,year,outdirname,m,mass))

if __name__ == "__main__":

    if not sys.warnoptions:
        import warnings
        warnings.simplefilter("ignore")
        os.environ["PYTHONWARNINGS"] = "ignore" # Also affect subprocesses

    #ex. python reloadFromPkl.py --indir cards/Mar05/01/ --pkl hadhadModel
    parser = argparse.ArgumentParser()
    parser.add_argument('--indir',      dest='indir',      default="indir",      help="indir")
    parser.add_argument('--year',       dest='year',       default="2017",       help="year")
    parser.add_argument('--outdir',     dest='outdir',     default="outdir",     help="outdir", nargs='+')
    parser.add_argument('--pkl',        dest='pkl',        default=["hadhadModel","hadmuModel","hadelModel"],    help="pickle name", nargs='+')
    parser.add_argument('--verbose',    dest='verbose',    action="store_true",  help="verbose")
    parser.add_argument('--addargs',    dest='addargs',    default="",  help='addargs')
    parser.add_argument('--mode',       dest='mode',       default="Significance",  help='mode')
    parser.add_argument('--expsig',     dest='expsig',     default="1",  help='expsig')
    parser.add_argument('--masses',     dest='masses',     default=['125'],  help='masses', nargs='+')
    parser.add_argument('--full',       dest='full',       action="store_true",  help='full')
    parser.add_argument('--dohtt',      dest='dohtt',      action="store_true",  help='dohtt')
    parser.add_argument('--mask',       dest='mask',       action="store_true",  help='mask')
    parser.add_argument('--nopass',     dest='nopass',     action="store_true",  help='nopass')
    args = parser.parse_args()

    reloadPkl(args)
