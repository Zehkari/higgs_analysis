import ROOT
from ROOT import TMath
import numpy as np
import matplotlib.pyplot as plt
%jsroot on

#Histogram properties
canvas = ROOT.TCanvas("Canvas","cz",800,600)
hist = ROOT.TH1F("h_M_Hyy","Diphoton invariant-mass ; m_{yy} [GeV] ; Events",60,100,160)
hist2 = ROOT.TH1F("bkg","Diphoton invariant-mass ; m_{yy} [GeV] ; Events",60,100,160)
#(variable name, title; x-axis label; y-axis label,min range,bin count,max range)

#Files
fA = "https://atlas-opendata.web.cern.ch/atlas-opendata/samples/2020/GamGam/Data/data_A.GamGam.root"
fB = "https://atlas-opendata.web.cern.ch/atlas-opendata/samples/2020/GamGam/Data/data_B.GamGam.root"
fC = "https://atlas-opendata.web.cern.ch/atlas-opendata/samples/2020/GamGam/Data/data_C.GamGam.root"
fD = "https://atlas-opendata.web.cern.ch/atlas-opendata/samples/2020/GamGam/Data/data_D.GamGam.root"
files = [fA,fB,fC,fD]

#Variables
Photon_1 = ROOT.TLorentzVector()
Photon_2 = ROOT.TLorentzVector()

#Data sortation
for name in files:
    f = ROOT.TFile.Open(name)
    tree = f.Get("mini")
    for event in tree:
        
        #Checking the trigger    
        if(tree.trigP):
            goodphoton_index = [0]*5
            goodphoton_n = 0
            photon_index = 0
            #Looping the photons per event
            for j in range(tree.photon_n):
                #Request of quality of the photon
                if(tree.photon_isTightID[j]):
                    #Transverse momentum and eta pre-selection
                    if(tree.photon_pt[j] > 25000 and (TMath.Abs(tree.photon_eta[j]) < 2.37)\
                       and (TMath.Abs(tree.photon_eta[j]) < 1.37 or TMath.Abs(tree.photon_eta[j]) > 1.52)):
                        goodphoton_n += 1  #count
                        goodphoton_index[photon_index]=j
                        photon_index += 1
            
            #Using the two selected photons
            if(goodphoton_n==2):
                goodphoton1_index = goodphoton_index[0]
                goodphoton2_index = goodphoton_index[1]
                #Isolation photon #1
                if((tree.photon_ptcone30[goodphoton1_index]/tree.photon_pt[goodphoton1_index] < 0.065)\
                   and (tree.photon_etcone20[goodphoton1_index] / tree.photon_pt[goodphoton1_index] < 0.065)):
                    #Isolation photon #2
                    if((tree.photon_ptcone30[goodphoton2_index]/tree.photon_pt[goodphoton2_index] < 0.065)\
                       and (tree.photon_etcone20[goodphoton2_index] / tree.photon_pt[goodphoton2_index] < 0.065)):
                        ##
                        Photon_1.SetPtEtaPhiE(tree.photon_pt[goodphoton1_index]/1000., tree.photon_eta[goodphoton1_index],\
                                          tree.photon_phi[goodphoton1_index],tree.photon_E[goodphoton1_index]/1000.)
                        Photon_2.SetPtEtaPhiE(tree.photon_pt[goodphoton2_index]/1000., tree.photon_eta[goodphoton2_index],\
                                          tree.photon_phi[goodphoton2_index],tree.photon_E[goodphoton2_index]/1000.)
                        
                        #Adding the two TLorentz vectors
                        Photon_12 = Photon_1 + Photon_2
                        #Filling histogram with the mass of the gamma-gamma system
                        hist.Fill(Photon_12.M())
                        
                        #Slicing data around a peak section
                        if Photon_12.M() < 115 or Photon_12.M() > 135:
                            hist2.Fill(Photon_12.M())

#Plotting data histogram
hist.Draw("E")
canvas.Draw()

#Plotting data without peak
hist2.Draw("E")
canvas.Draw()

#Fit bkg data with polynomial
hist2.Fit("pol3")
fit = hist2.GetFunction("pol3")
p0 = fit.GetParameter(0)
p1 = fit.GetParameter(1)
p2 = fit.GetParameter(2)
p3 = fit.GetParameter(3)
print(p0,p1,p2,p3)

#Sort data to python
E_index = []
data_index = []
no_bins = hist.GetNbinsX()
steps = 60/no_bins
bkg_index = []
for i in range(no_bins):
    E = 100+steps*i
    bkg = p0+(p1*E)+(p2*E**2)+(p3*E**3)
    bkg_index.append(bkg)
    data = hist.GetBinContent(i+1)
    E_index.append(E)
    data_index.append(data)

#fitting a line
plt.plot(E_index,data_index,".")
plt.plot(E_index,bkg_index,"r-")
plt.show()

#Sort data to python
E_index = []
data_index = []
no_bins = hist.GetNbinsX()
steps = 60/no_bins
bkg_index = []
for i in range(no_bins):
    E = 100+steps*i
    bkg = p0+(p1*E)+(p2*E**2)+(p3*E**3)
    bkg_index.append(bkg)
    data = hist.GetBinContent(i+1)
    E_index.append(E)
    data_index.append(data)

#plotting a line over raw data
plt.plot(E_index,data_index,".")
plt.plot(E_index,bkg_index,"r-")
plt.show()

#plot of peak
data_minus_bkg = []
for j in range(no_bins):
    data_minus_bkg.append(data_index[j]-bkg_index[j])

plt.plot(E_index,data_minus_bkg,".")
