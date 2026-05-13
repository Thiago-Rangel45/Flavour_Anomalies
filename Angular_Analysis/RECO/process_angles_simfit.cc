#include <iostream>
#include <vector>
#include <cmath>
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TString.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooWorkspace.h"
#include "RooArgSet.h"
#include "RooAbsData.h"

using namespace RooFit;

void process_angles_simfit() {
    RooAbsData::setDefaultStorageType(RooAbsData::Tree);

    TFile* f_in = TFile::Open("Filtered_RECO_Kinematics_NoFSR.root", "READ");
    if (!f_in || f_in->IsZombie()) {
        std::cerr << "Erro ao abrir Filtered_RECO_Kinematics_NoFSR.root" << std::endl;
        return;
    }
    TTree* tree = (TTree*)f_in->Get("Events");

    // Lendo como float (Correção anterior)
    float mass_Kpi, cos_theta_k, cos_theta_l, phi_angle, mll_fullfit;

    tree->SetBranchAddress("BToTrkTrkMuMu_fit_mass_Kpi", &mass_Kpi);
    tree->SetBranchAddress("BToTrkTrkMuMu_cos_theta_k", &cos_theta_k);
    tree->SetBranchAddress("BToTrkTrkMuMu_cos_theta_l", &cos_theta_l);
    tree->SetBranchAddress("BToTrkTrkMuMu_phi_angle", &phi_angle);
    tree->SetBranchAddress("BToTrkTrkMuMu_mll_fullfit", &mll_fullfit);

    // VARIÁVEIS EXATAMENTE COMO O SIMFIT ESPERA (INCLUINDO 'weight')
    RooRealVar ctK("ctK", "cos(#theta_{K})", -1, 1);
    RooRealVar ctL("ctL", "cos(#theta_{l})", -1, 1);
    RooRealVar phi("phi", "#phi", -3.14159, 3.14159);
    RooRealVar mass("mass", "m(#mu#muK#pi)", 5.0, 5.6, "GeV");
    RooRealVar rand("rand", "rand", 0, 1);
    RooRealVar weight("weight", "weight", 1.0); // Simfit exige essa variável
    RooArgSet argset(ctK, ctL, phi, rand, mass, weight);

    const int nBins = 9;
    std::vector<RooDataSet*> base_datasets;
    for (int i = 0; i < nBins; ++i) {
        base_datasets.push_back(new RooDataSet(Form("tmp_b%d", i), "tmp", argset));
    }

    TRandom3 rand_gen(12345);
    Long64_t nEntries = tree->GetEntries();
    std::cout << "Processando " << nEntries << " candidatos de B0 RECO..." << std::endl;

    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);

        double q2 = mll_fullfit * mll_fullfit;
        int b = -1;
        if (q2 >= 1.1 && q2 < 2.0) b = 0; 
        else if (q2 >= 2.0 && q2 < 4.3) b = 1;
        else if (q2 >= 4.3 && q2 < 6.0) b = 2; 
        else if (q2 >= 6.0 && q2 < 8.68) b = 3;
        else if (q2 >= 8.68 && q2 < 10.09) b = 4;
        else if (q2 >= 10.09 && q2 < 12.86) b = 5;
        else if (q2 >= 12.86 && q2 < 14.18) b = 6;
        else if (q2 >= 14.18 && q2 < 16.0) b = 7;
        else if (q2 >= 16.0 && q2 < 19.0) b = 8;
        
        if (b == -1) continue;

        if (std::isnan(cos_theta_k) || std::isnan(cos_theta_l) || std::isnan(phi_angle) || std::isnan(mass_Kpi)) continue;
        if (mass_Kpi < 5.0 || mass_Kpi > 5.6) continue;
        
        ctK.setVal(cos_theta_k); 
        ctL.setVal(cos_theta_l); 
        phi.setVal(phi_angle); 
        mass.setVal(mass_Kpi); 
        rand.setVal(rand_gen.Uniform(1));
        weight.setVal(1.0); // Peso sempre 1 para MC normal
        
        base_datasets[b]->add(argset);
    }
    f_in->Close();

    // Exportação - 1 ÚNICO ARQUIVO POR BIN COM OS 2 WORKSPACES!
    int year = 2016; // O ano padrão que o simfit vai procurar
    
    for (int i = 0; i < nBins; ++i) {
        if (base_datasets[i]->numEntries() == 0) continue;
        
        // Nome EXATO que o simfit vai procurar quando localFiles = true
        TString fileName = Form("./recoMCDataset_b%d_%d.root", i, year);
        TFile* f_out = TFile::Open(fileName, "RECREATE");

        RooDataSet* data_ev = (RooDataSet*)base_datasets[i]->reduce("rand < 0.5");
        RooDataSet* data_od = (RooDataSet*)base_datasets[i]->reduce("rand >= 0.5");

        // Salvar Paridade 0 (Even)
        RooWorkspace ws_ev(Form("ws_b%dp0", i));
        ws_ev.import(*data_ev, Rename(Form("dataset_b%dp0", i))); // Nome padrão utils.cc
        ws_ev.import(*data_ev, Rename("data")); // Backup name
        ws_ev.Write();

        // Salvar Paridade 1 (Odd)
        RooWorkspace ws_od(Form("ws_b%dp1", i));
        ws_od.import(*data_od, Rename(Form("dataset_b%dp1", i))); // Nome padrão utils.cc
        ws_od.import(*data_od, Rename("data")); // Backup name
        ws_od.Write();
        
        delete data_ev;
        delete data_od;
        
        f_out->Close();
        std::cout << "Bin " << i << " salvo com sucesso em " << fileName << std::endl;
    }
    
    for (int i = 0; i < nBins; ++i) { delete base_datasets[i]; }
}

int main() {
    process_angles_simfit();
    return 0;
}