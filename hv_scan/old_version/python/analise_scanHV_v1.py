
import ROOT             
from math import *

ndim = 1
nbin = 100
minE = [[0,100,40,200]]
#minE = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
maxE = [[100,400,80,400]]
#maxE = [[300,1000,200,600],[300,1000,200,600],[300,1000,200,600],[300,1000,200,600],[300,1000,200,600],[300,1000,200,600],[300,1000,200,600],[300,1000,200,600],[300,1000,200,600],[300,1000,200,600],[300,1000,200,600]]
    
def make_hists(name = " ",nbin=100,_minE0=0,_maxE0=1000,_minE1=0,_maxE1=1000,_minE2=0,_maxE2=1000,_minE3=0,_maxE3=1000) :
	h0 = ROOT.TH1D("h0","Energy", nbin, _minE0, _maxE0)
	h1 = ROOT.TH1D("h1","Energy", nbin, _minE1, _maxE1)
	h2 = ROOT.TH1D("h2","Energy", nbin, _minE2, _maxE2)
	h3 = ROOT.TH1D("h3","Energy", nbin, _minE3, _maxE3)
	f = ROOT.TFile(name)
	t = f.Get("T")
	iev = 0
	_hv0 = -1
	_hv1 = -1
	_hv2 = -1
	_hv3 = -1
	for event in t :
		if event.channel==0 : 
			h0.Fill(event.integral)
			h0.SetTitle("ch0 HV0")
		if event.channel==1 : 
			h1.Fill(event.integral)
			h1.SetTitle("ch1 HV1")
		if event.channel==2 : 
			h2.Fill(event.integral)
			h2.SetTitle("ch2 HV2")
		if event.channel==3 : 
			h3.Fill(event.integral)
			h3.SetTitle("ch3 HV3")
		if iev == 0 :
			_hv0 = event.hv4
			_hv1 = event.hv5
			_hv2 = event.hv6
			_hv3 = event.hv7
		iev += 1
	return (h0,_hv0,h1,_hv1,h2,_hv2,h3,_hv3)

hist0 = []
hist1 = []
hist2 = []
hist3 = []

namedata = [
            "HVCalibration/Grupo1-2/processed/mx_b_20170519_2044/mx_b_20170519_2044_000000.root"
            #"HVCalibration/Grupo1-2/processed/calibration/mx_b_20170519_2101/mx_b_20170519_2101_000000.root"
            #"HVCalibration/Grupo1-2/processed/mx_b_20170519_2119/mx_b_20170519_2119_000000.root",
            #"HVCalibration/Grupo1-2/processed/mx_b_20170519_2135/mx_b_20170519_2135_000000.root",
            #"HVCalibration/Grupo1-2/processed/mx_b_20170519_2149/mx_b_20170519_2149_000000.root",
            ##"HVCalibration/Grupo1-2/processed/mx_b_20170519_2219/mx_b_20170519_2219_000000.root",
            #"HVCalibration/Grupo1-2/processed/mx_b_20170523_1922/mx_b_20170523_1922_000000.root",
            #"HVCalibration/Grupo1-2/processed/mx_b_20170523_1958/mx_b_20170523_1958_000000.root",
            #"HVCalibration/Grupo1-2/processed/mx_b_20170523_2016/mx_b_20170523_2016_000000.root",
            #"HVCalibration/Grupo1-2/processed/mx_b_20170523_2042/mx_b_20170523_2042_000000.root",
            #"HVCalibration/Grupo1-2/processed/mx_b_20170523_2059/mx_b_20170523_2059_000000.root"
            ]

#canv = ROOT.TCanvas("canv","plots Energy", 1000, 1000 )
canvases = []
for i in range(0,ndim)	:
	canv = ROOT.TCanvas( ("canv_%d" %i),"plots Energy", 1000, 1000 ) 
	canvases.append( canv )
	h0,hv0,h1,hv1,h2,hv2,h3,hv3 = make_hists(namedata[i],nbin,minE[i][0],maxE[i][0],minE[i][1],maxE[i][1],minE[i][2],maxE[i][2],minE[i][3],maxE[i][3] ) #tuple
	#h0,hv0,h1,hv1,h2,hv2,h3,hv3 = make_hists(namedata[i]) #tuple
	hist0.append(h0)
	hist1.append(h1)
	hist2.append(h2)
	hist3.append(h3)
	print(hv0,hv1,hv2,hv3)

	hist0[i].SetLineColor(1)
	hist0[i].SetLineWidth(3)
	hist0[i].SetMarkerColor(9)
	hist0[i].SetFillColor(8)
	hist0[i].SetFillStyle(3001)
	hist0[i].SetMarkerStyle(20)
	hist0[i].GetXaxis().SetTitle("Energy")
	hist0[i].GetYaxis().SetTitle("events")

	hist1[i].SetLineColor(1)
	hist1[i].SetLineWidth(3)
	hist1[i].SetMarkerColor(2)
	hist1[i].SetFillColor(8)
	hist1[i].SetFillStyle(3001)
	hist1[i].SetMarkerStyle(20)
	hist1[i].GetXaxis().SetTitle("Energy")
	hist1[i].GetYaxis().SetTitle("events")

	hist2[i].SetLineColor(1)
	hist2[i].SetLineWidth(3)
	hist2[i].SetMarkerColor(2)
	hist2[i].SetFillColor(8)
	hist2[i].SetFillStyle(3001)
	hist2[i].SetMarkerStyle(20)
	hist2[i].GetXaxis().SetTitle("Energy")
	hist2[i].GetYaxis().SetTitle("events")

	hist3[i].SetLineColor(1)
	hist3[i].SetLineWidth(3)
	hist3[i].SetMarkerColor(2)
	hist3[i].SetFillColor(8)
	hist3[i].SetFillStyle(3001)
	hist3[i].SetMarkerStyle(20)
	hist3[i].GetXaxis().SetTitle("Energy")
	hist3[i].GetYaxis().SetTitle("events")

	#buf = "canvas %d".format(i)
	#canv = ROOT.TCanvas(buf,"plots Energy", 1000, 1000 )
	canvases[-1].Divide(2,2)

	canvases[-1].cd(1)
	hist0[i].Draw()
	canvases[-1].cd(2)
	hist1[i].Draw()
	canvases[-1].cd(3)
	hist2[i].Draw()
	canvases[-1].cd(4)
	hist3[i].Draw()

#canv.SaveAs("test.pdf")


