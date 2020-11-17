import ROOT
from ROOT import TMath
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
%jsroot on

#histogram properties
canvas = ROOT.TCanvas("Canvas","cz",800,600)
hist = ROOT.TH1F("h_M_Hyy","Diphoton invariant-mass ; m_{yy} [GeV] ; Events",30,100,160)
hist2 = ROOT.TH1F("bkg","Diphoton invariant-mass ; m_{yy} [GeV] ; Events",30,100,160)
#(variable name, title; x-axis label; y-axis label,bin count,min range,max range)

#files
fA = "https://atlas-opendata.web.cern.ch/atlas-opendata/samples/2020/GamGam/Data/data_A.GamGam.root"
fB = "https://atlas-opendata.web.cern.ch/atlas-opendata/samples/2020/GamGam/Data/data_B.GamGam.root"
fC = "https://atlas-opendata.web.cern.ch/atlas-opendata/samples/2020/GamGam/Data/data_C.GamGam.root"
fD = "https://atlas-opendata.web.cern.ch/atlas-opendata/samples/2020/GamGam/Data/data_D.GamGam.root"
files = [fA,fB,fC,fD]

#variables
Photon_1 = ROOT.TLorentzVector()
Photon_2 = ROOT.TLorentzVector()

#data sortation
for name in files:
    f = ROOT.TFile.Open(name)
    tree = f.Get("mini")
    for event in tree:
        
        #checking the trigger    
        if(tree.trigP):
            goodphoton_index = [0]*5
            goodphoton_n = 0
            photon_index = 0
            #looping the photons per event
            for j in range(tree.photon_n):
                #request of quality of the photon
                if(tree.photon_isTightID[j]):
                    #transverse momentum and eta pre-selection
                    if(tree.photon_pt[j] > 25000 and (TMath.Abs(tree.photon_eta[j]) < 2.37)\
                       and (TMath.Abs(tree.photon_eta[j]) < 1.37 or TMath.Abs(tree.photon_eta[j]) > 1.52)):
                        goodphoton_n += 1  #count
                        goodphoton_index[photon_index]=j
                        photon_index += 1
            
            #using the two selected photons
            if(goodphoton_n==2):
                goodphoton1_index = goodphoton_index[0]
                goodphoton2_index = goodphoton_index[1]
                #isolation photon #1
                if((tree.photon_ptcone30[goodphoton1_index]/tree.photon_pt[goodphoton1_index] < 0.065)\
                   and (tree.photon_etcone20[goodphoton1_index] / tree.photon_pt[goodphoton1_index] < 0.065)):
                    #isolation photon #2
                    if((tree.photon_ptcone30[goodphoton2_index]/tree.photon_pt[goodphoton2_index] < 0.065)\
                       and (tree.photon_etcone20[goodphoton2_index] / tree.photon_pt[goodphoton2_index] < 0.065)):
                        ##
                        Photon_1.SetPtEtaPhiE(tree.photon_pt[goodphoton1_index]/1000., tree.photon_eta[goodphoton1_index],\
                                          tree.photon_phi[goodphoton1_index],tree.photon_E[goodphoton1_index]/1000.)
                        Photon_2.SetPtEtaPhiE(tree.photon_pt[goodphoton2_index]/1000., tree.photon_eta[goodphoton2_index],\
                                          tree.photon_phi[goodphoton2_index],tree.photon_E[goodphoton2_index]/1000.)
                        
                        #adding the two TLorentz vectors
                        Photon_12 = Photon_1 + Photon_2
                        #filling histogram with the mass of the gamma-gamma system
                        hist.Fill(Photon_12.M())
                        
                        #slicing data around a peak section
                        if Photon_12.M() < 114 or Photon_12.M() > 134:
                            hist2.Fill(Photon_12.M())

#plotting data histogram
hist.Draw("E")
canvas.Draw()

#plotting data without peak
hist2.Draw("E")
canvas.Draw()

#fit bkg data with polynomial
hist2.Fit("pol6")
fit = hist2.GetFunction("pol6")
p0 = fit.GetParameter(0)
p1 = fit.GetParameter(1)
p2 = fit.GetParameter(2)
p3 = fit.GetParameter(3)
p4 = fit.GetParameter(4)
p5 = fit.GetParameter(5)
p6 = fit.GetParameter(6)
#p7 = fit.GetParameter(7)

#sort c++ data to python for analysis
E_index = []
data_index = []
no_bins = hist.GetNbinsX()
steps = 60/no_bins
bkg_index = []
for i in range(no_bins):
    E = 100+steps*i
    bkg = p0+(p1*E)+(p2*E**2)+(p3*E**3)+(p4*E**4)+(p5*E**5)+(p6*E**6)
    bkg_index.append(bkg)
    data = hist.GetBinContent(i+1)
    E_index.append(E)
    data_index.append(data)

#+(p4*E**4)+(p5*E**5)+(p6*E**6)

#fitting a line
plt.plot(E_index,data_index,".",label='Data')
plt.plot(E_index,bkg_index,"r-",label = 'Background Fit')
plt.xlabel('mγγ [GeV]')
plt.ylabel('Events[Bins]')
plt.rcParams["figure.figsize"] = (8,2)
plt.legend()
plt.show()

#subtracting peak from bkg
data_minus_bkg = []
for j in range(no_bins):
    data_minus_bkg.append(data_index[j]-bkg_index[j])

#gaussian
def gauss(x, a, mean, sigma):
    return a * np.exp(-(x - mean) ** 2 / (2 * sigma ** 2))
	
#fit gaussian
popt,pcov = curve_fit(gauss,E_index,data_minus_bkg,p0=[150,125,2])
a, mean, sigma = popt
print(a,mean,sigma)

#plotting gaussian and data
xp = np.linspace(100,160,100)
plt.plot(E_index,data_minus_bkg,".",label = 'Data - Background')
plt.plot(xp,gauss(xp,a,mean,sigma),'r--', label = 'Gaussian Fit')
plt.xlabel("mγγ [GeV]")
plt.ylabel("Events[Bins]")
plt.rcParams["figure.figsize"] = (10,6)
plt.legend()
plt.show()

#finding local maximum
imax = np.argmax(gauss(xp,a,mean,sigma))
print("Local maxima y-cord: {}".format(imax))
print("Local maxima x-cord: {}".format(xp[imax]))

#invariant mass peak
print("The inavariant mass value for the peak is: {} GeV".format(round(xp[imax],1)))