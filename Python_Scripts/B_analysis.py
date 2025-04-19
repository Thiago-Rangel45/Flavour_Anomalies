import ROOT
from ROOT import RooRealVar, RooDataHist, RooPlot, RooCrystalBall, RooArgList
import numpy as np



INPUT_DIR = "/afs/cern.ch/user/t/tdeandra/Flavour_Analysis/data/"   
INPUT_FILE = "reco_ntuple_LMNR_1.root"
TREE_NAME = "ntuple"

MU_MASS = 0.10566
PION_MASS = 0.140
KAON_MASS = 0.494

def create_histogram(name, title, bins, x_min, x_max, x_label, y_label="N_{Events}"):
    hist = ROOT.TH1F(name, title, bins, x_min, x_max)
    hist.GetXaxis().SetTitle(x_label)
    hist.GetYaxis().SetTitle(y_label)
    return hist

def create_lorentz(pt, eta, phi, mass):
    vec = ROOT.TLorentzVector()
    vec.SetPtEtaPhiM(pt, eta, phi, mass)
    return vec

def draw_and_save(hist, canvas_name, filename, color=ROOT.kRed, marker=20):
    canvas = ROOT.TCanvas(canvas_name, "", 800, 600)
    hist.SetMarkerStyle(marker)
    hist.SetLineColor(color)
    hist.Draw("HIST")
    canvas.SaveAs(filename)

def normalize(hist):
    integral = hist.Integral()
    if integral > 0:
        hist.Scale(1.0 / integral)

file = ROOT.TFile.Open(INPUT_DIR + INPUT_FILE)
tree = file.Get(TREE_NAME)

h_kstar  = create_histogram("h_kstmass", "", 10, 0, 10, "mass (K^{*0} #rightarrow K^{+} #pi^{-}) (GeV)")
h_kstbar = create_histogram("h_kstbarmass", "", 10, 0, 10, "mass (K^{*0}_{bar} #rightarrow K^{-} #pi^{+}) (GeV)")
h_dimu   = create_histogram("h_dimu", "", 10, 0, 10, "mass (#mu^{+} #mu^{-}) (GeV)")
h_B      = create_histogram("h_Bmass", "", 20, 4.9, 5.7, "mass (B^{0} #rightarrow K^{*0} #mu^{+} #mu^{-}) (GeV)")
h_Bbar   = create_histogram("h_Bbarmass", "", 20, 4.9, 5.7, "mass (B^{0}_{bar} #rightarrow K^{*0}_{bar} #mu^{+} #mu^{-}) (GeV)")


branches = {
    "mumPt": np.zeros(1), "mumEta": np.zeros(1), "mumPhi": np.zeros(1),
    "mupPt": np.zeros(1), "mupEta": np.zeros(1), "mupPhi": np.zeros(1),
    "kstTrkmPt": np.zeros(1), "kstTrkmEta": np.zeros(1), "kstTrkmPhi": np.zeros(1),
    "kstTrkpPt": np.zeros(1), "kstTrkpEta": np.zeros(1), "kstTrkpPhi": np.zeros(1),
}

for branch, arr in branches.items():
    tree.SetBranchAddress(branch, arr)


for i in range(tree.GetEntries()):
    tree.GetEntry(i)

    # Reconstrução dos vetores de momento
    mu1 = create_lorentz(branches["mumPt"][0], branches["mumEta"][0], branches["mumPhi"][0], MU_MASS)
    mu2 = create_lorentz(branches["mupPt"][0], branches["mupEta"][0], branches["mupPhi"][0], MU_MASS)
    dimu = mu1 + mu2
    h_dimu.Fill(dimu.M())

    # K*0 (K+ pi-)
    pion_m = create_lorentz(branches["kstTrkmPt"][0], branches["kstTrkmEta"][0], branches["kstTrkmPhi"][0], PION_MASS)
    kaon_p = create_lorentz(branches["kstTrkpPt"][0], branches["kstTrkpEta"][0], branches["kstTrkpPhi"][0], KAON_MASS)
    kstar = pion_m + kaon_p
    h_kstar.Fill(kstar.M())

    # K*0_bar (K- pi+)
    pion_p = create_lorentz(branches["kstTrkpPt"][0], branches["kstTrkpEta"][0], branches["kstTrkpPhi"][0], PION_MASS)
    kaon_m = create_lorentz(branches["kstTrkmPt"][0], branches["kstTrkmEta"][0], branches["kstTrkmPhi"][0], KAON_MASS)
    kstbar = pion_p + kaon_m
    h_kstbar.Fill(kstbar.M())

    # Reconstrução do B ou B_bar
    if abs(kstar.M() - 0.896) < abs(kstbar.M() - 0.896):
        if 2.9 < dimu.M() < 3.3 and 0.9 < kstar.M() < 1.2:
            B = kstar + dimu
            h_B.Fill(B.M())
    else:
        if 2.9 < dimu.M() < 3.3 and 0.9 < kstbar.M() < 1.2:
            B_bar = kstbar + dimu
            h_Bbar.Fill(B_bar.M())

for hist in [h_kstar, h_kstbar, h_dimu, h_Bbar]:
    normalize(hist)

draw_and_save(h_kstar,  "c1", "Kstar_mass.png")
draw_and_save(h_kstbar, "c2", "Kstar_bar_mass.png")
draw_and_save(h_dimu,   "c3", "dimuon_mass.png")
draw_and_save(h_B,      "c4", "B_mass.png")
draw_and_save(h_Bbar,   "c5", "Bbar_mass.png")