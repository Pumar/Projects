
import ROOT             
from math import *
    
nbin = 100
minE = 0
maxE = 300

hintegral1 = ROOT.TH1D("hintegral1","Energy", nbin, minE, maxE) 
hintegral2 = ROOT.TH1D("hintegral2","Energy", nbin, minE, maxE) 
hint = [] 

cut      = ["channel==0","channel==1"]
title    = ["HV = 100","HV = 200"]
namedata = ["HVCalibration/Group1/files_Cs137_group1_1000/mx_b_20170407_1736/mx_b_20170407_1736_000000.root",
            "HVCalibration/Group1/files_Cs137_group1_1000/mx_b_20170407_1736/mx_b_20170407_1736_000000.root"]

canv = ROOT.TCanvas("canv","plots Energy", 1000,1000 )
canv.Divide(1,2) 

#for i in range(2) :
	
f = ROOT.TFile(namedata[0])
t = f.Get("T")
t.Project("hintegral1","integral",cut[0])
hintegral1.SetTitle(title[0]) 
hintegral1.SetLineColor(2)
canv.cd(1)
hintegral1.Draw()
f.Close()

f = ROOT.TFile(namedata[1])
t = f.Get("T")
t.Project("hintegral2","integral",cut[1])
hintegral2.SetTitle(title[1]) 
hintegral2.SetLineColor(2)
canv.cd(2)
hintegral2.Draw()

canv.SaveAs("plots_test.pdf") 

