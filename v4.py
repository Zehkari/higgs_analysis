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
fA = "data_A.GamGam.root"
fB = "data_B.GamGam.root"
fC = "data_C.GamGam.root"
fD = "data_D.GamGam.root"
files = [fA,fB,fC,fD]

#variables
Photon_1 = ROOT.TLorentzVector()
Photon_2 = ROOT.TLorentzVector()

#data selection
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

#plotting data histogram
hist.Draw("E")
canvas.Draw()

#sort c++ data to python for analysis
E_index = []
data_index = []
sideband_E = []
sideband_data = []
no_bins = hist.GetNbinsX()
steps = 60/no_bins
for i in range(no_bins):
    E = 100+steps*i
    data = hist.GetBinContent(i+1)
    E_index.append(E)
    data_index.append(data)
    if (E < 115 or E > 135):
        sideband_E.append(E)
        sideband_data.append(data)

#data errors
data_ers = []
for i in data_index:
    N = np.sqrt(i)
    data_ers.append(N)

print(data_index)
print("\n\n",data_ers)

#sideband errors
side_ers = []
for i in sideband_data:
    N = np.sqrt(i)
    side_ers.append(N)

#data errors
data_ers = []
for i in data_index:
    N = np.sqrt(i)
    data_ers.append(N)

#sideband errors
side_ers = []
for i in sideband_data:
    N = np.sqrt(i)
    side_ers.append(N)

#sliced data plot
plt.plot(sideband_E,sideband_data,'.')

#polynomial function
def pol(x,p0,p1,p2,p3):
    return p0+(p1*x)+(p2*x**2)+(p3*x**3)

#fit polynomial
popt,pcov = curve_fit(pol,sideband_E,sideband_data,sigma=side_ers,absolute_sigma=True)
p0,p1,p2,p3 = popt
print("p0 = ",p0,"\np1 = ",p1,"\np2 = ",p2,"\np3 = ",p3)
print("\n",pcov)

#data minus bkg
bkg_index = pol(np.array(E_index),p0,p1,p2,p3)
data_minus_bkg = np.array(data_index-bkg_index)

#polynomial errors (apprroximation)
#diag of covariance matrix
perr = np.sqrt(np.sum(np.diag(pcov)))

#background errors
bkg_ers = []
for i in data_ers:
    bkg_ers.append(perr)

#plotting background fit
plt.plot(E_index,data_index,".",label='Data')
plt.errorbar(E_index,data_index, yerr=data_ers, fmt='.',
                    ecolor='lightgray', elinewidth=2, capsize=0,label='Errors')
plt.plot(E_index,bkg_index,"r-",label = 'Background Fit')
plt.xlabel('mγγ [GeV]')
plt.ylabel('Events[Bins]')
plt.legend()
plt.show()

#data minus bkg errors
#propagation of error
minus_ers = ((np.array(data_ers))**2 + (np.array(bkg_ers))**2)**(1/2)

#gaussian function
def gauss(x, a, mean, sigma):
    return a * np.exp(-(x - mean) ** 2 / (2 * sigma ** 2))

#fit gaussian
max_y = max(data_minus_bkg)
popt_g,pcov_g = curve_fit(gauss,E_index,data_minus_bkg,p0=[max_y,125,2])
a, mean, sigma = popt_g
print("amp value: ",a,"\nMean: ",mean,"\nStd: ",sigma)

#plotting gaussian and data
xp = np.linspace(100,160,100)
plt.plot(E_index,data_minus_bkg,".",label = 'Data')
plt.errorbar(E_index,data_minus_bkg, yerr=minus_ers, fmt='.',
                    ecolor='lightgray', elinewidth=1, capsize=4,label='Errors')
plt.plot(xp,gauss(xp,a,mean,sigma),'r-', label = 'Gaussian Fit')
plt.xlabel("mγγ [GeV]")
plt.ylabel("Data - bkg")
plt.legend(loc='upper right')
plt.show()
plt.savefig("gauss.png",dpi=300)

#local maxima y-cord
print("The local maxima y-cord: {}".format(round(a,1)))

#invariant mass peak
print("The inavariant mass value for the peak is: {} GeV".format(round(mean,1)))