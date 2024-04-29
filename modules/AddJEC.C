#include <vector>
#include <string>

#include <TString.h>
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"

void AddJEC(std::string inputFilePath, std::string outputFilePath, std::string treeName, std::string jecVersion){
  std::cout << "=========================================================" << std::endl;
  std::cout << "AddJEC()::Start" << std::endl;
  std::cout << "inputFilePath: " << inputFilePath << std::endl;
  std::cout << "outputFilePath: " << outputFilePath << std::endl;
  std::cout << "treeName: " << treeName << std::endl;
  //
  //
  //
  std::cout << "Using following JEC files:" << std::endl;
  std::string inDir("CondFormats/JetMETObjects/data/");
  std::string txtFileRes(inDir+jecVersion+"_L2L3Residual_AK4PFPuppi.txt");
  std::string txtFileL3(inDir+jecVersion+"_L3Absolute_AK4PFPuppi.txt");
  std::string txtFileL2(inDir+jecVersion+"_L2Relative_AK4PFPuppi.txt");
  std::string txtFileL1(inDir+jecVersion+"_L1FastJet_AK4PFPuppi.txt");
  std::cout << txtFileRes << std::endl;
  std::cout << txtFileL3 << std::endl;
  std::cout << txtFileL2 << std::endl;
  std::cout << txtFileL1 << std::endl;

  std::cout << "Setup JetCorrectorParameters " << std::endl;

  JetCorrectorParameters ResJetPar(txtFileRes);
  JetCorrectorParameters L3JetPar(txtFileL3);
  JetCorrectorParameters L2JetPar(txtFileL2);
  JetCorrectorParameters L1JetPar(txtFileL1);

  //
  // Load the JetCorrectorParameter objects into a vector, IMPORTANT: THE ORDER MATTERS HERE !!!!
  //
  std::cout << "Setup std::vector<JetCorrectorParameters> " << std::endl;
  std::vector<JetCorrectorParameters> vPar;
  vPar.push_back(ResJetPar);
  vPar.push_back(L3JetPar);
  vPar.push_back(L2JetPar);
  vPar.push_back(L1JetPar);

  std::cout << "Setup FactorizedJetCorrector " << std::endl;
  FactorizedJetCorrector jetCorrector(vPar);

  std::cout << "Setup std::vector<JetCorrectorParameters> for L1Only " << std::endl;
  std::vector<JetCorrectorParameters> vParL1;
  vParL1.push_back(L1JetPar);

  std::cout << "Setup FactorizedJetCorrector " << std::endl;
  FactorizedJetCorrector jetCorrectorL1(vParL1);

  //
  //
  //
  std::cout << "Setup inputFile " << std::endl;
  TFile inputTreeFile(inputFilePath.c_str(), "OPEN");

  std::cout << "Retrieve inputTree " << std::endl;
  // TTree* inTree = (TTree*) inputTreeFile.Get(treeName.c_str());
  TTree* inTree = inputTreeFile.Get<TTree>(treeName.c_str());

  if (!inTree) {
    std::cout << "No TTree. Skip this file!" << std::endl;
    return;
  }

  auto nEntries = inTree->GetEntries();
  std::cout << "Retrieve input tree which has nentries: " << nEntries << std::endl;

  // if (nEntries == 0): return;
  //==================================================
  //
  //
  //
  //==================================================
  TFile outputTreeFile(outputFilePath.c_str(), "RECREATE");

  inTree->SetBranchStatus("Jet_pt",0);
  inTree->SetBranchStatus("Jet_mass",0);

  TTree* outTree = inTree->CloneTree(0);

  const Int_t nJetsMax = 200;
  Float_t Jet_pt_nano[nJetsMax];
  Float_t Jet_pt_raw[nJetsMax];
  Float_t Jet_pt_updated[nJetsMax];

  Float_t Jet_mass_nano[nJetsMax];
  Float_t Jet_mass_raw[nJetsMax];
  Float_t Jet_mass_updated[nJetsMax];

  Float_t Jet_rawFactor_nano[nJetsMax];
  Float_t Jet_rawFactor_updated[nJetsMax];

  Float_t PuppiMet_pt_updated;
  Float_t PuppiMet_phi_updated;

  Float_t PuppiMet_pt_nano;
  Float_t PuppiMet_phi_nano;
  Float_t PuppiMet_sumEt_nano;

  outTree->Branch("Jet_pt_nano",        Jet_pt_nano,           "Jet_pt_nano[nJet]/F");
  outTree->Branch("Jet_mass_nano",      Jet_mass_nano,         "Jet_mass_nano[nJet]/F");
  outTree->Branch("Jet_pt_raw",         Jet_pt_raw,            "Jet_pt_raw[nJet]/F");
  outTree->Branch("Jet_mass_raw",       Jet_mass_raw,          "Jet_mass_raw[nJet]/F");
  outTree->Branch("Jet_rawFactor_nano", Jet_rawFactor_nano,    "Jet_rawFactor_nano[nJet]/F");
  outTree->Branch("Jet_rawFactor",      Jet_rawFactor_updated, "Jet_rawFactor[nJet]/F");

  // TBranch* b_jet_pt = outTree->GetBranch("Jet_pt");
  // TBranch* b_jet_mass = outTree->GetBranch("Jet_mass");
  // b_jet_pt->SetAddress(Jet_pt_updated);
  // b_jet_mass->SetAddress(Jet_mass_updated);
  outTree->Branch("Jet_pt",             Jet_pt_updated,        "Jet_pt[nJet]/F");
  outTree->Branch("Jet_mass",           Jet_mass_updated,      "Jet_mass[nJet]/F");

  outTree->Branch("PuppiMet_pt",         &PuppiMet_pt_updated,    "PuppiMet_pt/F");
  outTree->Branch("PuppiMet_phi",        &PuppiMet_phi_updated,   "PuppiMet_phi/F");

  outTree->Branch("PuppiMet_pt_nano",    &PuppiMet_pt_nano,    "PuppiMet_pt_nano/F");
  outTree->Branch("PuppiMet_phi_nano",   &PuppiMet_phi_nano,   "PuppiMet_phi_nano/F");
  outTree->Branch("PuppiMet_sumEt_nano", &PuppiMet_sumEt_nano, "PuppiMet_sumEt_nano/F");


  //
  //
  //
  //
  inTree->SetBranchStatus("Jet_pt",1);
  inTree->SetBranchStatus("Jet_mass",1);

  TTreeReader reader(inTree);

  TTreeReaderValue<Float_t> Rho_fixedGridRhoFastjetAll = {reader, "Rho_fixedGridRhoFastjetAll"};
  TTreeReaderValue<Int_t> nJet = {reader, "nJet"};
  TTreeReaderArray<Float_t> Jet_pt   = {reader, "Jet_pt"};
  TTreeReaderArray<Float_t> Jet_phi  = {reader, "Jet_phi"};
  TTreeReaderArray<Float_t> Jet_eta  = {reader, "Jet_eta"};
  TTreeReaderArray<Float_t> Jet_mass = {reader, "Jet_mass"};
  TTreeReaderArray<Float_t> Jet_rawFactor = {reader, "Jet_rawFactor"};
  TTreeReaderArray<Float_t> Jet_area = {reader, "Jet_area"};
  TTreeReaderArray<Float_t> Jet_chEmEF = {reader, "Jet_chEmEF"};
  TTreeReaderArray<Float_t> Jet_neEmEF = {reader, "Jet_neEmEF"};
  TTreeReaderArray<Float_t> Jet_muonSubtrFactor = {reader, "Jet_muonSubtrFactor"};

  TTreeReaderValue<Int_t>   nCorrT1METJet       = {reader, "nCorrT1METJet"};
  TTreeReaderArray<Float_t> CorrT1METJet_rawPt  = {reader, "CorrT1METJet_rawPt"};
  TTreeReaderArray<Float_t> CorrT1METJet_eta    = {reader, "CorrT1METJet_eta"};
  TTreeReaderArray<Float_t> CorrT1METJet_phi    = {reader, "CorrT1METJet_phi"};
  TTreeReaderArray<Float_t> CorrT1METJet_area   = {reader, "CorrT1METJet_area"};
  TTreeReaderArray<Float_t> CorrT1METJet_EmEF   = {reader, "CorrT1METJet_EmEF"};
  TTreeReaderArray<Float_t> CorrT1METJet_muonSubtrFactor = {reader, "CorrT1METJet_muonSubtrFactor"};

  TTreeReaderValue<Float_t> RawPuppiMET_pt    = {reader, "RawPuppiMET_pt"};
  TTreeReaderValue<Float_t> RawPuppiMET_phi   = {reader, "RawPuppiMET_phi"};

  TTreeReaderValue<Float_t> PuppiMET_pt    = {reader, "PuppiMET_pt"};
  TTreeReaderValue<Float_t> PuppiMET_phi   = {reader, "PuppiMET_phi"};
  TTreeReaderValue<Float_t> PuppiMET_sumEt = {reader, "PuppiMET_sumEt"};

  for(uint iEntry=0; iEntry < nEntries; iEntry++){
    reader.SetLocalEntry(iEntry);
    inTree->GetEntry(iEntry);

    if (iEntry % 1000 == 0){
      std::cout << iEntry << " / " << nEntries << std::endl;
    }

    float rho = *Rho_fixedGridRhoFastjetAll;

    float met_px_updated = (*RawPuppiMET_pt) * TMath::Cos((*RawPuppiMET_phi));
    float met_py_updated = (*RawPuppiMET_pt) * TMath::Sin((*RawPuppiMET_phi));

    //==================================================
    //
    //
    //
    //==================================================
    for (int iJet=0; iJet < *nJet; iJet++){
      float recojet_pt        = Jet_pt[iJet];
      float recojet_eta       = Jet_eta[iJet];
      float recojet_phi       = Jet_phi[iJet];
      float recojet_mass      = Jet_mass[iJet];
      float recojet_area      = Jet_area[iJet];
      float recojet_rawFactor = Jet_rawFactor[iJet];
      float recojet_muonSubtrFactor = Jet_muonSubtrFactor[iJet];
      float recojet_pt_raw    = (1.f - recojet_rawFactor) * recojet_pt;
      float recojet_mass_raw  = (1.f - recojet_rawFactor) * recojet_mass;

      //
      //
      //
      jetCorrector.setJetPt(recojet_pt_raw);
      jetCorrector.setJetEta(recojet_eta);
      jetCorrector.setJetPhi(recojet_phi);
      jetCorrector.setJetA(recojet_area);
      jetCorrector.setRho(rho);
      float jec = jetCorrector.getCorrection();
      float recojet_pt_updated   = recojet_pt_raw * jec;
      float recojet_mass_updated = recojet_mass_raw * jec;

      Jet_pt_nano[iJet] = recojet_pt;
      Jet_mass_nano[iJet] = recojet_mass;
      Jet_rawFactor_nano[iJet] = recojet_rawFactor;

      Jet_pt_raw[iJet]    = recojet_pt_raw;
      Jet_mass_raw[iJet]  = recojet_mass_raw;

      Jet_pt_updated[iJet]   = recojet_pt_updated;
      Jet_mass_updated[iJet] = recojet_mass_updated;
      Jet_rawFactor_updated[iJet] = 1.f - (recojet_pt_raw/recojet_pt_updated);

      //========================================================================
      //
      // Type-I MET corrections
      //
      //========================================================================
      float recojet_cosphi = TMath::Cos(recojet_phi);
      float recojet_sinphi = TMath::Sin(recojet_phi);

      //
      // 1. Get the jet raw pt without muons included
      //
      float recojet_pt_noMuRaw = recojet_pt_raw * (1.f - recojet_muonSubtrFactor);

      //
      // 2. Apply JEC on this muon-less jet raw pt
      //
      jetCorrector.setJetPt(recojet_pt_noMuRaw);
      jetCorrector.setJetEta(recojet_eta);
      jetCorrector.setJetPhi(recojet_phi);
      jetCorrector.setJetA(recojet_area);
      jetCorrector.setRho(rho);
      float jecForNoMuRaw = jetCorrector.getCorrection();
      float recojet_pt_noMuL1L2L3 = recojet_pt_raw * jecForNoMuRaw;
      float recojet_px_noMuL1L2L3 = recojet_pt_noMuL1L2L3 * recojet_cosphi;
      float recojet_py_noMuL1L2L3 = recojet_pt_noMuL1L2L3 * recojet_sinphi;

      //
      // 3. 
      //
      float recojet_pt_noMuOnlyL1 = recojet_pt_noMuRaw;
      jetCorrectorL1.setJetPt(recojet_pt_raw);
      jetCorrectorL1.setJetEta(recojet_eta);
      jetCorrectorL1.setJetPhi(recojet_phi);
      jetCorrectorL1.setJetA(recojet_area);
      jetCorrectorL1.setRho(rho);
      float jecOnlyL1ForNoMuRaw = jetCorrectorL1.getCorrection();
      recojet_pt_noMuOnlyL1 *= jecOnlyL1ForNoMuRaw;
      float recojet_px_noMuOnlyL1 = recojet_pt_noMuOnlyL1 * TMath::Cos(recojet_phi);
      float recojet_py_noMuOnlyL1 = recojet_pt_noMuOnlyL1 * TMath::Sin(recojet_phi);

      if ((recojet_pt_noMuL1L2L3 > 15.f) && (Jet_chEmEF[iJet] + Jet_neEmEF[iJet] < 0.9)){
        met_px_updated -= (recojet_px_noMuL1L2L3 - recojet_px_noMuOnlyL1);
        met_py_updated -= (recojet_py_noMuL1L2L3 - recojet_py_noMuOnlyL1);
      }
    }

    //
    //
    //
    for (int iJet=0; iJet < *nCorrT1METJet; iJet++){
      float recojet_pt_raw    = CorrT1METJet_rawPt[iJet];
      float recojet_eta       = CorrT1METJet_eta[iJet];
      float recojet_phi       = CorrT1METJet_phi[iJet];
      float recojet_area      = CorrT1METJet_area[iJet];
      float recojet_muonSubtrFactor = CorrT1METJet_muonSubtrFactor[iJet];
      //========================================================================
      //
      // Type-I MET corrections
      //
      //========================================================================
      float recojet_cosphi = TMath::Cos(recojet_phi);
      float recojet_sinphi = TMath::Sin(recojet_phi);

      //
      // 1. Get the jet raw pt without muons included
      //
      float recojet_pt_noMuRaw = recojet_pt_raw * (1.f - recojet_muonSubtrFactor);

      //
      // 2. Apply JEC on this muon-less jet raw pt
      //
      jetCorrector.setJetPt(recojet_pt_noMuRaw);
      jetCorrector.setJetEta(recojet_eta);
      jetCorrector.setJetPhi(recojet_phi);
      jetCorrector.setJetA(recojet_area);
      jetCorrector.setRho(rho);
      float jecForNoMuRaw = jetCorrector.getCorrection();
      float recojet_pt_noMuL1L2L3 = recojet_pt_raw * jecForNoMuRaw;
      float recojet_px_noMuL1L2L3 = recojet_pt_noMuL1L2L3 * recojet_cosphi;
      float recojet_py_noMuL1L2L3 = recojet_pt_noMuL1L2L3 * recojet_sinphi;

      //
      // 3.
      //
      float recojet_pt_noMuOnlyL1 = recojet_pt_noMuRaw;
      jetCorrectorL1.setJetPt(recojet_pt_raw);
      jetCorrectorL1.setJetEta(recojet_eta);
      jetCorrectorL1.setJetPhi(recojet_phi);
      jetCorrectorL1.setJetA(recojet_area);
      jetCorrectorL1.setRho(rho);
      float jecOnlyL1ForNoMuRaw = jetCorrectorL1.getCorrection();
      recojet_pt_noMuOnlyL1 *= jecOnlyL1ForNoMuRaw;
      float recojet_px_noMuOnlyL1 = recojet_pt_noMuOnlyL1 * TMath::Cos(recojet_phi);
      float recojet_py_noMuOnlyL1 = recojet_pt_noMuOnlyL1 * TMath::Sin(recojet_phi);

      if ((recojet_pt_noMuL1L2L3 > 15.f) && (CorrT1METJet_EmEF[iJet] < 0.9)){
        met_px_updated -= (recojet_px_noMuL1L2L3 - recojet_px_noMuOnlyL1);
        met_py_updated -= (recojet_py_noMuL1L2L3 - recojet_py_noMuOnlyL1);
      }
    }

    PuppiMet_pt_updated  = TMath::Sqrt((met_px_updated * met_px_updated) + (met_py_updated * met_py_updated));
    PuppiMet_phi_updated = TMath::ATan2(met_px_updated, met_py_updated);

    PuppiMet_pt_nano    = *(PuppiMET_pt);
    PuppiMet_phi_nano   = *(PuppiMET_phi);
    PuppiMet_sumEt_nano = *(PuppiMET_sumEt);
    outTree->Fill();
  }

  std::cout << "NEntries in input  :" << nEntries << std::endl;
  std::cout << "NEntries in output :" << outTree->GetEntries() << std::endl;

  outTree->Write();
  outputTreeFile.Close();
  inputTreeFile.Close();
  std::cout << "AddJEC()::Finished" << std::endl;
}
