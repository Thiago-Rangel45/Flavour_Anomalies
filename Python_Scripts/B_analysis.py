import ROOT
import numpy as np

input_dir = "/afs/cern.ch/user/t/tdeandra/Flavour_Analysis/data/"   
input_file = "ntuple_flat_LMNR_PostRefitMomenta_test_2022F_skimSoftMu_1.root"
tree_name = "ntuple"

root_file = ROOT.TFile.Open(input_dir + input_file)
tree = root_file.Get(tree_name)

h1 = ROOT.TH1F("h_kstmass", "", 10, 0, 10) # Histograma para o K*
h1.GetXaxis().SetTitle("mass (K^{*0} -> K^{+} #pi^{-}) (GeV)")
h1.GetYaxis().SetTitle("N_{Events}")

h2 = ROOT.TH1F("h_kstbarmass", "", 10, 0, 10)     # Histograma para o K* bar
h2.GetXaxis().SetTitle("mass (K^{*0}_{bar} -> K^{-} #pi^{+}) (GeV)")
h2.GetYaxis().SetTitle("N_{Events}")

h3 = ROOT.TH1F("h_dimu", "", 10, 0, 10)
h3.GetXaxis().SetTitle("mass {#mu^{+} #mu^{-}} GeV")   # Histograma para o dimu
h3.GetYaxis().SetTitle("N_{Events}")

h4 = ROOT.TH1F("h_Bmass", "", 80, 4.9, 5.7)
h4.GetXaxis().SetTitle("mass (B^{0} -> K^{0*} #mu^{+} #mu^{-}) (GeV)")   # Histograma para o B
h4.GetYaxis().SetTitle("N_{Events}")

h5 = ROOT.TH1F("h_Bbarmass", "", 80, 4.9, 5.7)
h5.GetXaxis().SetTitle("mass (B^{0}_{bar} -> K^{0*}_{bar} #mu^{+} #mu^{-}) (GeV)")   # Histograma para o B bar
h5.GetYaxis().SetTitle("N_{Events}")

# Variáveis dos múons e K* e K* Tracker + e -
mu_mass = 0.10566       # massa de repouso do múon em GeV
pion_mass = 0.140       # massa de repouso do píon em GeV
kaon_mass = 0.494       # massa de repouso do káon em GeV

mupEta = np.zeros(1, dtype=float)
mupPhi = np.zeros(1, dtype=float)
mupPt = np.zeros(1, dtype=float)

mumEta = np.zeros(1, dtype=float)
mumPhi = np.zeros(1, dtype=float)
mumPt = np.zeros(1, dtype=float)

kstTrkmPt = np.zeros(1, dtype=float)
kstTrkmEta = np.zeros(1, dtype=float)
kstTrkmPhi = np.zeros(1, dtype=float)

kstTrkpPt = np.zeros(1, dtype=float)
kstTrkpEta = np.zeros(1, dtype=float)
kstTrkpPhi = np.zeros(1, dtype=float)

# SetBranchAddress
tree.SetBranchAddress("mumPt", mumPt)
tree.SetBranchAddress("mumEta", mumEta)
tree.SetBranchAddress("mumPhi", mumPhi)

tree.SetBranchAddress("mupPt", mupPt)
tree.SetBranchAddress("mupEta", mupEta)
tree.SetBranchAddress("mupPhi", mupPhi)

tree.SetBranchAddress("kstTrkmPt", kstTrkmPt)
tree.SetBranchAddress("kstTrkmEta", kstTrkmEta)
tree.SetBranchAddress("kstTrkmPhi", kstTrkmPhi)

tree.SetBranchAddress("kstTrkpPt", kstTrkpPt)
tree.SetBranchAddress("kstTrkpEta", kstTrkpEta)
tree.SetBranchAddress("kstTrkpPhi", kstTrkpPhi)

mu1 = ROOT.TLorentzVector() # mu -
mu2 = ROOT.TLorentzVector() # mu +
dimu = ROOT.TLorentzVector() # dimu

pion_m = ROOT.TLorentzVector() # Pion + 
kaon_p = ROOT.TLorentzVector() # K -
kstar = ROOT.TLorentzVector()   # K*

pion_p = ROOT.TLorentzVector() # Pion - 
kaon_m = ROOT.TLorentzVector() # K +
kstar_bar = ROOT.TLorentzVector() # K* bar

B_cand = ROOT.TLorentzVector() # B
B_bar_cand = ROOT.TLorentzVector() # B bar

n_entries = tree.GetEntries()
for i in range(tree.GetEntries()):
    tree.GetEntry(i)
    
    # seleção dos K* e K* bar
    pion_m.SetPtEtaPhiM(kstTrkmPt[0], kstTrkmEta[0], kstTrkmPhi[0], pion_mass) # Lorentz Vector para o Pion -
    kaon_p.SetPtEtaPhiM(kstTrkpPt[0], kstTrkpEta[0], kstTrkpPhi[0], kaon_mass) # Lorentz Vector para o Kaon +
    kstar = pion_m + kaon_p
    h1.Fill(kstar.M())

    pion_p.SetPtEtaPhiM(kstTrkpPt[0], kstTrkpEta[0], kstTrkpPhi[0], pion_mass) # Pion +
    kaon_m.SetPtEtaPhiM(kstTrkmPt[0], kstTrkmEta[0], kstTrkmPhi[0], kaon_mass) # Kaon -
    kstar_bar = pion_p + kaon_m
    h2.Fill(kstar_bar.M())

    # seleção dos múons
    mu1.SetPtEtaPhiM(mumPt[0], mumEta[0], mumPhi[0], mu_mass)
    mu2.SetPtEtaPhiM(mupPt[0], mupEta[0], mupPhi[0], mu_mass)
    dimu = mu1 + mu2
    h3.Fill(dimu.M())

    if abs(kstar.M() - 0.896) < abs(kstar_bar.M() - 0.896): # seleção para que seja um B0
        if 0.9 < kstar.M() < 1.2 and 2.9 < dimu.M() < 3.3:  # selecionando somente na janela de massa do K* e J/psi
            B_cand = kstar + dimu
            h4.Fill(B_cand.M())
    else:
        if 0.9 < kstar_bar.M() < 1.2 and 2.9 < dimu.M() < 3.3:  # selecionando somente na janela de massa do K* bar e J/psi
            B_bar_cand = kstar_bar + dimu
            h5.Fill(B_bar_cand.M())

# Desenho do histograma da massa do K
c1 = ROOT.TCanvas("c1", "", 800, 600)
h1.SetMarkerStyle(22)
h1.Draw("HIST")  
c1.SaveAs("K_mass_distribution.png")

# Desenho do histograma de massa do K*
c2 = ROOT.TCanvas("c2", "", 800, 600)
h2.SetMarkerStyle(20)
h2.SetLineColor(ROOT.kRed)
h2.Draw("HIST")
c2.SaveAs("K_bar_mass_window.png")

# Desenho do histograma de massa do dimu
c3= ROOT.TCanvas("c3", "", 800, 600)
h3.SetMarkerStyle(20)
h3.SetLineColor(ROOT.kRed)
h3.Draw("HIST")
c3.SaveAs("dimu_mass_window.png")

c4 = ROOT.TCanvas("c4", "", 800, 600)
h4.SetMarkerStyle(20)
h4.SetLineColor(ROOT.kRed)
h4.Draw("HIST")
c4.SaveAs("B_mass_window.png")

c5 = ROOT.TCanvas("c5", "", 800, 600)
h5.SetMarkerStyle(20)
h5.SetLineColor(ROOT.kRed)
h5.Draw("HIST")
c5.SaveAs("B_bar_mass_window.png")