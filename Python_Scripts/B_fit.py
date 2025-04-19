import ROOT
import numpy as np
from ROOT import RooRealVar, RooDataHist, RooPlot, RooCrystalBall, RooArgList

# Constantes
MU_MASS   = 0.105658  # Massa do múon em GeV
PION_MASS = 0.139570  # Massa do píon em GeV
KAON_MASS = 0.493677  # Massa do kaon em GeV

# Diretório e arquivo ROOT
INPUT_DIR = "/afs/cern.ch/user/t/tdeandra/Flavour_Analysis/data/"
INPUT_FILE = "reco_ntuple_LMNR_1.root"
TREE_NAME = "ntuple"

# Abrindo o arquivo e a árvore
file = ROOT.TFile.Open(INPUT_DIR + INPUT_FILE)
tree = file.Get(TREE_NAME)

# Inicializando o histograma
hist = ROOT.TH1F("hist", "B^{0} mass", 20, 4.9, 5.7)
hist.GetXaxis().SetTitle("mass (K^{*0} #rightarrow K^{+} #pi^{-}) (GeV)")
hist.GetYaxis().SetTitle("N_{Events}")
hist.SetMarkerStyle(20)
hist.SetLineColor(ROOT.kBlue)

# Definindo os branches
branches = {
    "mumPt": np.zeros(1), "mumEta": np.zeros(1), "mumPhi": np.zeros(1),
    "mupPt": np.zeros(1), "mupEta": np.zeros(1), "mupPhi": np.zeros(1),
    "kstTrkmPt": np.zeros(1), "kstTrkmEta": np.zeros(1), "kstTrkmPhi": np.zeros(1),
    "kstTrkpPt": np.zeros(1), "kstTrkpEta": np.zeros(1), "kstTrkpPhi": np.zeros(1),
}

# Ligando os branches
for branch, arr in branches.items():
    tree.SetBranchAddress(branch, arr)

# Função para preencher o histograma
def fill_histogram(tree, hist):
    # Inicializando os quatro-vetores
    mu1 = ROOT.TLorentzVector()
    mu2 = ROOT.TLorentzVector()
    pion_m = ROOT.TLorentzVector()
    kaon_p = ROOT.TLorentzVector()
    pion_p = ROOT.TLorentzVector()
    kaon_m = ROOT.TLorentzVector()

    # Loop para preencher o histograma
    for i in range(tree.GetEntries()):
        tree.GetEntry(i)

        # Definindo as massas de cada partícula
        mu1.SetPtEtaPhiM(branches["mumPt"][0], branches["mumEta"][0], branches["mumPhi"][0], MU_MASS)
        mu2.SetPtEtaPhiM(branches["mupPt"][0], branches["mupEta"][0], branches["mupPhi"][0], MU_MASS)
        dimu = mu1 + mu2

        pion_m.SetPtEtaPhiM(branches["kstTrkmPt"][0], branches["kstTrkmEta"][0], branches["kstTrkmPhi"][0], PION_MASS)
        kaon_p.SetPtEtaPhiM(branches["kstTrkpPt"][0], branches["kstTrkpEta"][0], branches["kstTrkpPhi"][0], KAON_MASS)
        kstar = pion_m + kaon_p

        pion_p.SetPtEtaPhiM(branches["kstTrkpPt"][0], branches["kstTrkpEta"][0], branches["kstTrkpPhi"][0], PION_MASS)
        kaon_m.SetPtEtaPhiM(branches["kstTrkmPt"][0], branches["kstTrkmEta"][0], branches["kstTrkmPhi"][0], KAON_MASS)
        kstbar = pion_p + kaon_m

        # Filtrando os eventos conforme os critérios
        if abs(kstar.M() - 0.896) < abs(kstbar.M() - 0.896):
            if 2.9 < dimu.M() < 3.3 and 0.9 < kstar.M() < 1.2:
                B = kstar + dimu
                hist.Fill(B.M())

# Preenchendo o histograma
fill_histogram(tree, hist)

# Definindo as variáveis de RooFit
mass = RooRealVar("mass", "mass", 4.9, 5.7)
data_B = RooDataHist("data_B", "B^{0} mass data", RooArgList(mass), hist)

# PDFs
mean  = RooRealVar("mean",  "Mean", 5.28, 5.20, 5.35)
sigma = RooRealVar("sigma", "Sigma", 0.2, 0., 1.)
alpha = RooRealVar("alpha", "Alpha", 1.5, 0.5, 5.0)
n     = RooRealVar("n",     "n",     3.0, 0.5, 10.0)
cb_pdf = RooCrystalBall("cb_pdf", "Crystal Ball PDF", mass, mean, sigma, alpha, n)

p0 = RooRealVar("p0", "p0", 0.0, -1.0, 1.0)
background = ROOT.RooPolynomial("background", "Linear Background", mass, RooArgList(p0))

fsig = RooRealVar("fsig", "signal fraction", 0.9, 0.0, 1.0)
model = ROOT.RooAddPdf("model", "CB + Linear", RooArgList(cb_pdf, background), RooArgList(fsig))

# Ajustando o modelo aos dados
model.fitTo(data_B)
fit_result = model.fitTo(data_B, ROOT.RooFit.Save())

# Criando o gráfico
frame = mass.frame()
data_B.plotOn(frame, ROOT.RooFit.Name("data"))
model.plotOn(frame, ROOT.RooFit.Name("model"))
model.plotOn(frame, ROOT.RooFit.Components("background"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name("background"))
model.plotOn(frame, ROOT.RooFit.Components("cb_pdf"), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name("cb_pdf"))

frame.SetTitle("")
frame.GetXaxis().SetTitle("m(K^{+} #pi^{-} #mu^{+} #mu^{-}) (GeV)")
frame.GetYaxis().SetTitle("N_{Events} / (0.04 GeV)")

# Configuração do Canvas
c_fit = ROOT.TCanvas("c_fit", "Crystal Ball + Background Fit", 700, 700)
c_fit.SetTickx(1)
c_fit.SetTicky(1)
frame.Draw()

# Legenda
legend = ROOT.TLegend(0.12, 0.65, 0.43, 0.88)  
legend.SetBorderSize(0)
legend.SetFillStyle(0)
legend.SetTextFont(42)
legend.SetTextSize(0.035)

legend.AddEntry(frame.findObject("data"), "Data", "lep")
legend.AddEntry(frame.findObject("model"), "Fit (CB + linear)", "l")
legend.AddEntry(frame.findObject("cb_pdf"), "Crystal Ball", "l")
legend.AddEntry(frame.findObject("background"), "Background", "l")
legend.Draw()

# Adicionando as anotações no gráfico
latex = ROOT.TLatex()
latex.SetNDC()
latex.SetTextFont(62)
latex.SetTextSize(0.045)
latex.DrawLatex(0.10, 0.91, "CMS Simulation")

latex.SetTextFont(42)
latex.SetTextSize(0.030)
latex.DrawLatex(0.60, 0.85, "8.41 < q^{2} < 10.89 GeV^{2}")
latex.DrawLatex(0.60, 0.80, "0.9 < K^{0*} < 1.2 GeV^{2}")
latex.DrawLatex(0.60, 0.75, "#chi^{2}/ndf = 11.13")

# Salvando a imagem
c_fit.SaveAs("B_mass_fit.png")

n_param = fit_result.floatParsFinal().getSize()
chi2 = frame.chiSquare(n_param)
print(chi2)
print(n_param)