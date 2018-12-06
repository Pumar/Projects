Histo_training_S.SetLineColor(2)
Histo_training_S.SetMarkerColor(2)
Histo_training_S.SetFillColor(2)
Histo_testing_S.SetLineColor(2)
Histo_testing_S.SetMarkerColor(2)
Histo_testing_S.SetFillColor(2)
 
Histo_training_B.SetLineColor(4)
Histo_training_B.SetMarkerColor(4)
Histo_training_B.SetFillColor(4)
Histo_testing_B.SetLineColor(4)
Histo_testing_B.SetMarkerColor(4)
Histo_testing_B.SetFillColor(4)
 
# Histogram fill styles
Histo_training_S.SetFillStyle(3001)
Histo_training_B.SetFillStyle(3001)
Histo_testing_S.SetFillStyle(0)
Histo_testing_B.SetFillStyle(0)
 
# Histogram marker styles
Histo_testing_S.SetMarkerStyle(20)
Histo_testing_B.SetMarkerStyle(20)
 
# Set titles
Histo_training_S.GetXaxis().SetTitle("Classifier, SVM [rbf kernel, C=1, gamma=0.005]")
Histo_training_S.GetYaxis().SetTitle("Counts/Bin")
 
# Draw the objects
c1 = ROOT.TCanvas("c1","",800,600)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
Histo_training_S.Draw("HIST")
Histo_training_B.Draw("HISTSAME")
Histo_testing_S.Draw("EPSAME")
Histo_testing_B.Draw("EPSAME")
 
# Reset the y-max of the plot
ymax = max([h.GetMaximum() for h in [Histo_training_S,Histo_training_B,Histo_testing_S,Histo_testing_B] ])
ymax *=1.2
Histo_training_S.SetMaximum(ymax)
 
# Create Legend
c1.cd(1).BuildLegend( 0.42,  0.72,  0.57,  0.88).SetFillColor(0)
 
# Add custom title
l1=ROOT.TLatex()
l1.SetNDC();
l1.DrawLatex(0.26,0.93,"Classification with TMVA (ROOT)")
 
