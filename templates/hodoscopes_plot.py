import ROOT
from ctypes import *
import argparse
import glob

ROOT.gSystem.Load("lib/libH4Analysis.so")
ROOT.gROOT.SetBatch(True)

parser = argparse.ArgumentParser (description = 'make plots from hodoscopes 1 and 2')
parser.add_argument('-r', '--run' , help='run to be processed', type = int, default = 15153)
parser.add_argument('--path' , default = '/eos/user/c/cbasile/ECAL_TB2021/ntuples_templates/ntuples_v1')
parser.add_argument('--outdir' , default = '/eos/user/c/cbasile/ECAL_TB2021/ntuples_templates/ntuples_v1')
parser.add_argument('--crystal' , default = 'C2')
parser.add_argument('--tag' ,   )
parser.add_argument('--energy' , default = 100, type = int)
parser.add_argument('--hcx' , default = 0., type = float)
parser.add_argument('--hcy' , default = 0., type = float)
parser.add_argument('--tcx' , default = 0., type = float)
parser.add_argument('--tcy' , default = 0., type = float)

args = parser.parse_args ()
run = args.run
path = args.path
outdir = args.outdir 

track_plot = False

h4 = ROOT.TChain('h4')
nfiles = h4.Add(f'{path}/{run}/*.root')
print(f' + Added {nfiles} files')

ADCmax = args.energy*32 # 32 ADC/GeV

h4.Draw('amp_max['+args.crystal+']:h1Y.clusters.Y_:h1X.clusters.X_>>hodo1(40, -20, 20, 40, -20, 20, 0, 10000)', ' n_clusters>=1', 'PROFCOLZ')
h_hodo1 = ROOT.gPad.GetPrimitive("hodo1")
h1_area = ROOT.TBox(args.hcx-2, args.hcy-2, args.hcx+2, args.hcy+2)
h4.Draw('amp_max[%s]:h1X.clusters.X_>>hodo1_x(40,-20,20, 1000,0, %d)'%(args.crystal, ADCmax), 'n_clusters>=1', 'COLZ')
h_hodo1_x = ROOT.gPad.GetPrimitive("hodo1_x")
h4.Draw('amp_max[%s]:h1Y.clusters.Y_>>hodo1_y(40,-20,20, 500,0, %d)'%(args.crystal, ADCmax), 'n_clusters>=1', 'COLZ')
h_hodo1_y = ROOT.gPad.GetPrimitive("hodo1_y")

h4.Draw('amp_max['+args.crystal+']:h2Y.clusters.Y_:h2X.clusters.X_>>hodo2(40, -20, 20, 40, -20, 20,0, 10000)', ' n_clusters >= 1', 'PROFCOLZ')
h_hodo2 = ROOT.gPad.GetPrimitive("hodo2")

if track_plot :
   h4.Draw('amp_max['+args.crystal+']:Y:X>>tracks(40, -20, 20, 40, -20, 20,0, 10000)', ' n_clusters>=1', 'PROFCOLZ')
   h_tracks= ROOT.gPad.GetPrimitive("tracks")
   xmean = args.tcx #projX.GetBinCenter(projX.GetMaximumBin())
   ymean = args.tcy #projY.GetBinCenter(projY.GetMaximumBin())
   print(f' = track <X> = {xmean} mm <Y> = {ymean} mm')
   grmax= ROOT.TH1F('gr_%s_%s'%(args.crystal,args.energy),'gr_%s_%s'%(args.crystal,args.energy),1,xmean-1,xmean+1)
   grmax.SetBinContent(1,ymean)
   trk_area = ROOT.TBox(xmean-2, ymean-2, xmean+2, ymean+2)

c = ROOT.TCanvas("c", "", 1600, 800)
ROOT.gStyle.SetOptStat(0)

margin = 0.14
c.Divide(2,1)
c.cd(1)
ROOT.gPad.SetMargin(margin, margin, margin, margin)
h_hodo1.GetXaxis().SetTitle("X (mm)")
h_hodo1.GetYaxis().SetTitle("Y (mm)")
h_hodo1.SetTitle("HODOSCOPE 1")
h_hodo1.Draw("PROFCOLZ0")
h1_area.SetLineColor(ROOT.kRed)
h1_area.SetLineWidth(2)
h1_area.SetFillStyle(0)
h1_area.Draw("same")
c.cd(2)
ROOT.gPad.SetMargin(margin, margin, margin, margin)
h_hodo2.GetXaxis().SetTitle("X (mm)")
h_hodo2.GetYaxis().SetTitle("Y (mm)")
h_hodo2.SetTitle("HODOSCOPE 2")
h_hodo2.Draw("PROFCOLZ0")
c.SaveAs(outdir+'/AmpHodo_%s_%dGeV%s.png'%(args.crystal,args.energy,args.tag if(args.tag) else ''))

c_proj = ROOT.TCanvas("c_proj", "", 1600, 800)
c_proj.Divide(2,1)
c_proj.cd(1)
ROOT.gPad.SetMargin(margin, margin, margin, margin)
h_hodo1_x.GetXaxis().SetTitle("X (mm)")
h_hodo1_x.Draw("PROFCOLZ0")
c_proj.cd(2)
ROOT.gPad.SetMargin(margin, margin, margin, margin)
h_hodo1_y.Draw("COLZ0")
h_hodo1_y.GetXaxis().SetTitle("Y (mm)")
c_proj.SaveAs(outdir+'/AmpHodoX_%s_%dGeV%s.png'%(args.crystal,args.energy,args.tag if(args.tag) else ''))

if track_plot:
   c2 = ROOT.TCanvas("c2", "", 1024, 1024)
   ROOT.gPad.SetMargin(margin, margin, margin, margin)
   h_tracks.GetXaxis().SetTitle("X (mm)")
   h_tracks.GetYaxis().SetTitle("Y (mm)")
   h_tracks.SetTitle("TRACKS")
   h_tracks.Draw("PROFCOLZ0")
   grmax.SetMarkerColor(ROOT.kRed)
   grmax.SetMarkerStyle(34)
   grmax.SetMarkerSize(3)
   grmax.Draw("LPsame")
   trk_area.SetLineColor(ROOT.kRed)
   trk_area.SetFillStyle(0)
   trk_area.Draw("same")
   c2.SaveAs(outdir+'/AmpTracks_%s_run%s_%dGeV.png'%(args.crystal, run,args.energy)) 
