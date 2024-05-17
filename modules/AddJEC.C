#include <vector>
#include <string>
#include <numeric> // std::iota

#include <TString.h>
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TH1D.h>
#include <TH2D.h>

#include <Math/PtEtaPhiE4D.h>
#include <Math/PtEtaPhiM4D.h>
#include <Math/LorentzVector.h>

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"

FactorizedJetCorrector* SetupJetCorrector(std::string inDirPath, std::string L1="", std::string L2="", std::string L2L3Res=""){
  std::cout << "Setup FactorizedJetCorrector with the following JEC files:" << std::endl;
  std::vector<JetCorrectorParameters> vPar;
  std::string fullPath;
  if (L1 != ""){
    fullPath = Form("%s/%s.txt",inDirPath.c_str(), L1.c_str());
    std::cout << "Loading JEC:" << fullPath << std::endl;
    JetCorrectorParameters L1JetPar(fullPath);
    vPar.push_back(L1JetPar);
  }
  if (L2 != ""){
    fullPath = Form("%s/%s.txt",inDirPath.c_str(), L2.c_str());
    std::cout << "Loading JEC:" << fullPath << std::endl;
    JetCorrectorParameters L2JetPar(fullPath);
    vPar.push_back(L2JetPar);
  }
  if (L2L3Res != ""){
    fullPath = Form("%s/%s.txt",inDirPath.c_str(), L2L3Res.c_str());
    std::cout << "Loading JEC:" << fullPath << std::endl;
    JetCorrectorParameters ResJetPar(fullPath);
    vPar.push_back(ResJetPar);
  }

  FactorizedJetCorrector* jc = new FactorizedJetCorrector(vPar);
  return jc;
}

TH2D* SetupJetVetoMap(std::string inDirPath,std::string inFileName, std::string vetomapName){
  std::string fullPath = Form("%s/%s.root",inDirPath.c_str(), inFileName.c_str());

  TFile* inFile = nullptr;
  inFile = new TFile(fullPath.c_str(),"OPEN");

  TH2D* h2 = nullptr;
  if (inFile)
    h2 = (TH2D*)inFile->Get(vetomapName.c_str());

  TH2D* h2C = nullptr;
  if (h2){
    h2C = (TH2D*)h2->Clone((vetomapName+"_C").c_str());
    h2C->SetDirectory(0);
  }

  if (inFile)
    inFile->Close();

  return h2C;
}

void AddJEC(std::string inputFilePath, std::string outputFilePath, std::string treeName, std::string jecVersion, bool isMC=false){
  std::cout << "=========================================================" << std::endl;
  std::cout << "AddJEC()::Start" << std::endl;
  std::cout << "inputFilePath: " << inputFilePath << std::endl;
  std::cout << "outputFilePath: " << outputFilePath << std::endl;
  std::cout << "treeName: " << treeName << std::endl;

  //==================================================
  //
  // Setup
  //
  //==================================================
  std::string inDirPath("");
  FactorizedJetCorrector* jetCorrector = nullptr;
  TH2D* h2_jetvetomap = nullptr;

  // inDirPath("CondFormats/JetMETObjects/data/");
  // jetCorrector = SetupJetCorrector(inDirPath,
  //   jecVersion+"_L1FastJet_AK4PFPuppi",
  //   jecVersion+"_L2Relative_AK4PFPuppi",
  //   jecVersion+"_L2L3Residual_AK4PFPuppi"
  // );

  if (jecVersion == "Prompt24_V2M"){
    inDirPath = "CondFormats/JetMETObjects/data/Prompt24_V2M/";
    if(isMC){
      jetCorrector = SetupJetCorrector(inDirPath,
        "",
        "Winter24Run3_V1_MC_L2Relative_AK4PUPPI",
        ""
      );
    }else{
      jetCorrector = SetupJetCorrector(inDirPath,
        "",
        "Winter24Run3_V1_MC_L2Relative_AK4PUPPI",
        "Prompt24_Run2024BC_V2M_DATA_L2L3Residual_AK4PFPuppi"
      );
    }
    h2_jetvetomap = SetupJetVetoMap(
      inDirPath,
      "jetveto2024BC_V2M","jetvetomap"
    );
  }
  //==================================================
  //
  //
  //
  //==================================================
  bool applyL1OnNoMuP4 = false;
  FactorizedJetCorrector* jetCorrectorL1 = nullptr;
  if (applyL1OnNoMuP4){
    jetCorrectorL1 = SetupJetCorrector(inDirPath,
      jecVersion+"_L1FastJet_AK4PFPuppi","",""
    );
  }

  //==================================================
  //
  //
  //
  //==================================================
  std::cout << "Setup inputFile " << std::endl;
  TFile inputTreeFile(inputFilePath.c_str(), "OPEN");

  std::cout << "Retrieve inputTree " << std::endl;
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
  // Make output TTree
  //
  //==================================================
  TFile outputTreeFile(outputFilePath.c_str(), "RECREATE");

  // First, lets switch off these branches
  // in the input tree so that when we clone
  // the tree, the new tree doesn't have them
  // by default.
  inTree->SetBranchStatus("Jet_pt",0);
  inTree->SetBranchStatus("Jet_mass",0);
  inTree->SetBranchStatus("Jet_rawFactor",0);
  inTree->SetBranchStatus("PuppiMET_pt",0);
  inTree->SetBranchStatus("PuppiMET_phi",0);

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

  Int_t Jet_InVetoRegion[nJetsMax];

  Int_t Jet_IdxSortByPt[nJetsMax];

  Float_t Jet_pt_noMuL1L2L3[nJetsMax];
  Float_t CorrT1METJet_pt[nJetsMax];
  Float_t CorrT1METJet_pt_noMuL1L2L3[nJetsMax];

  Float_t PuppiMet_pt_updated;
  Float_t PuppiMet_phi_updated;

  Float_t PuppiMet_pt_nano;
  Float_t PuppiMet_phi_nano;
  Float_t PuppiMet_sumEt_nano;

  //
  // This is where we manually put in the branches
  // that we omitted when cloning original input tree
  //
  outTree->Branch("Jet_pt",              Jet_pt_updated,        "Jet_pt[nJet]/F");
  outTree->Branch("Jet_mass",            Jet_mass_updated,      "Jet_mass[nJet]/F");
  outTree->Branch("PuppiMET_pt",         &PuppiMet_pt_updated,  "PuppiMET_pt/F");
  outTree->Branch("PuppiMET_phi",        &PuppiMet_phi_updated, "PuppiMET_phi/F");
  // Add these extra "Jet" branches
  outTree->Branch("Jet_pt_nano",         Jet_pt_nano,           "Jet_pt_nano[nJet]/F");
  outTree->Branch("Jet_mass_nano",       Jet_mass_nano,         "Jet_mass_nano[nJet]/F");
  outTree->Branch("Jet_pt_raw",          Jet_pt_raw,            "Jet_pt_raw[nJet]/F");
  outTree->Branch("Jet_mass_raw",        Jet_mass_raw,          "Jet_mass_raw[nJet]/F");
  outTree->Branch("Jet_rawFactor_nano",  Jet_rawFactor_nano,    "Jet_rawFactor_nano[nJet]/F");
  outTree->Branch("Jet_rawFactor",       Jet_rawFactor_updated, "Jet_rawFactor[nJet]/F");
  outTree->Branch("Jet_InVetoRegion",    Jet_InVetoRegion,      "Jet_InVetoRegion[nJet]/I");

  outTree->Branch("Jet_IdxSortByPt",     Jet_IdxSortByPt,       "Jet_IdxSortByPt[nJet]/I");

  outTree->Branch("Jet_pt_noMuL1L2L3",   Jet_pt_noMuL1L2L3, "Jet_pt_noMuL1L2L3[nJet]/F");
  outTree->Branch("CorrT1METJet_pt",     CorrT1METJet_pt, "CorrT1METJet_pt[nCorrT1METJet]/F");
  outTree->Branch("CorrT1METJet_pt_noMuL1L2L3", CorrT1METJet_pt_noMuL1L2L3, "CorrT1METJet_pt_noMuL1L2L3[nCorrT1METJet]/F");

  outTree->Branch("PuppiMET_pt_nano",    &PuppiMet_pt_nano,     "PuppiMET_pt_nano/F");
  outTree->Branch("PuppiMET_phi_nano",   &PuppiMet_phi_nano,    "PuppiMET_phi_nano/F");
  outTree->Branch("PuppiMET_sumEt_nano", &PuppiMet_sumEt_nano,  "PuppiMET_sumEt_nano/F");

  //==================================================
  //
  // Setup TTreeReader for the input tree
  //
  //==================================================
  // We need to switch back on these branches in the input tree
  // so that we can retrieve the values stored in the tree.
  inTree->SetBranchStatus("Jet_pt",1);
  inTree->SetBranchStatus("Jet_mass",1);
  inTree->SetBranchStatus("Jet_rawFactor",1);
  inTree->SetBranchStatus("PuppiMET_pt",1);
  inTree->SetBranchStatus("PuppiMET_phi",1);

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

  TTreeReaderValue<Int_t>   nCorrT1METJet        = {reader, "nCorrT1METJet"};
  TTreeReaderArray<Float_t> CorrT1METJet_rawPt   = {reader, "CorrT1METJet_rawPt"};
  TTreeReaderArray<Float_t> CorrT1METJet_eta     = {reader, "CorrT1METJet_eta"};
  TTreeReaderArray<Float_t> CorrT1METJet_phi     = {reader, "CorrT1METJet_phi"};
  TTreeReaderArray<Float_t> CorrT1METJet_area    = {reader, "CorrT1METJet_area"};
  TTreeReaderArray<Float_t> CorrT1METJet_muonSubtrFactor = {reader, "CorrT1METJet_muonSubtrFactor"};
  TTreeReaderArray<Float_t> CorrT1METJet_EmEF    = {reader, "CorrT1METJet_EmEF"};
  bool hasEmEF = true;

  TTreeReaderValue<Float_t> RawPuppiMET_pt    = {reader, "RawPuppiMET_pt"};
  TTreeReaderValue<Float_t> RawPuppiMET_phi   = {reader, "RawPuppiMET_phi"};

  TTreeReaderValue<Float_t> PuppiMET_pt    = {reader, "PuppiMET_pt"};
  TTreeReaderValue<Float_t> PuppiMET_phi   = {reader, "PuppiMET_phi"};
  TTreeReaderValue<Float_t> PuppiMET_sumEt = {reader, "PuppiMET_sumEt"};

  bool firstEntry=true;

  for(uint iEntry=0; iEntry < nEntries; iEntry++){
    reader.SetLocalEntry(iEntry);
    inTree->GetEntry(iEntry);

    if (firstEntry){
      if (CorrT1METJet_EmEF.GetSetupStatus() < 0){
        hasEmEF = false;
        std::cout << "========================================================================" << std::endl;
        std::cout << "WARNING. CorrT1METJet_EmEF does not exist!!! Will not try to read branch" << std::endl;
        std::cout << "========================================================================" << std::endl;
      }
      firstEntry=false;
    }

    if (iEntry % 1000 == 0){
      std::cout << iEntry << " / " << nEntries << std::endl;
    }

    float rho = *Rho_fixedGridRhoFastjetAll;

    //==================================================
    //
    // For MET Type-1 re-correction with latest JEC.
    // Start of with the RawPuppiMET and we'll correct
    // each term jet-by-jet in the jet loop.
    //
    //==================================================
    float met_px_updated = (*RawPuppiMET_pt) * TMath::Cos((*RawPuppiMET_phi));
    float met_py_updated = (*RawPuppiMET_pt) * TMath::Sin((*RawPuppiMET_phi));

    //==================================================
    //
    // Loop over main jet collection
    //
    //==================================================
    //
    std::vector<float> v_jet_pt_updated(*nJet,1.f);

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
      // Get JEC.
      //
      jetCorrector->setJetPt(recojet_pt_raw);
      jetCorrector->setJetEta(recojet_eta);
      jetCorrector->setJetPhi(recojet_phi);
      jetCorrector->setJetA(recojet_area);
      jetCorrector->setRho(rho);
      float jec = jetCorrector->getCorrection();
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

      v_jet_pt_updated[iJet] = recojet_pt_updated;

      //========================================================================
      //
      // Check if jet is in jetveto regioon
      //
      //========================================================================
      Jet_InVetoRegion[iJet] = 0;
      if(h2_jetvetomap){
        int i1 = h2_jetvetomap->GetXaxis()->FindBin(recojet_eta);
        int j1 = h2_jetvetomap->GetYaxis()->FindBin(recojet_phi);
        Jet_InVetoRegion[iJet] = int(h2_jetvetomap->GetBinContent(i1, j1) > 0);
      }

      //========================================================================
      //
      // MET Type-1 correction
      //
      //========================================================================
      //
      // 1. Get the jet raw pt without muons included
      //
      float recojet_pt_noMuRaw   = recojet_pt_raw * (1.f - recojet_muonSubtrFactor);
      // This is not exactly the correct thing to do. We shouldn't be using the full
      // jet eta and phi but the eta and phi of muon-less jet p4. Not possible with
      // current NanoAODs.
      float recojet_eta_noMuRaw  = recojet_eta;
      float recojet_phi_noMuRaw  = recojet_phi;

      //
      // 2. Apply JEC on this muon-less jet raw p4
      // The JEC used is the same JEC we apply on the full raw p4
      float recojet_pt_noMuL1L2L3 = jec * recojet_pt_noMuRaw;
      float recojet_px_noMuL1L2L3 = recojet_pt_noMuL1L2L3 * TMath::Cos(recojet_phi_noMuRaw);
      float recojet_py_noMuL1L2L3 = recojet_pt_noMuL1L2L3 * TMath::Sin(recojet_phi_noMuRaw);
      Jet_pt_noMuL1L2L3[iJet] = recojet_pt_noMuL1L2L3;

      //
      // 3.
      //
      float recojet_pt_noMuOnlyL1 = recojet_pt_noMuRaw;
      float recojet_px_noMuOnlyL1 = recojet_pt_noMuRaw * TMath::Cos(recojet_phi_noMuRaw);
      float recojet_py_noMuOnlyL1 = recojet_pt_noMuRaw * TMath::Sin(recojet_phi_noMuRaw);
      if (applyL1OnNoMuP4 && jetCorrectorL1){
        jetCorrectorL1->setJetPt(recojet_pt_raw);
        jetCorrectorL1->setJetEta(recojet_eta);
        jetCorrectorL1->setJetPhi(recojet_phi);
        jetCorrectorL1->setJetA(recojet_area);
        jetCorrectorL1->setRho(rho);
        float jecOnlyL1ForNoMuRaw = jetCorrectorL1->getCorrection();
        recojet_pt_noMuOnlyL1 *= jecOnlyL1ForNoMuRaw;
        recojet_px_noMuOnlyL1 = recojet_pt_noMuL1L2L3 * TMath::Cos(recojet_phi_noMuRaw);
        recojet_py_noMuOnlyL1 = recojet_pt_noMuL1L2L3 * TMath::Sin(recojet_phi_noMuRaw);
      }
      float recojet_emEF = Jet_chEmEF[iJet] + Jet_neEmEF[iJet];

      //
      // 4. Calculate the Type-1 correction and add it to the Raw MET
      //
      if ((recojet_pt_noMuL1L2L3 > 15.f) && (recojet_emEF < 0.9f)){
        met_px_updated -= (recojet_px_noMuL1L2L3 - recojet_px_noMuOnlyL1);
        met_py_updated -= (recojet_py_noMuL1L2L3 - recojet_py_noMuOnlyL1);
      }
    }

    //==================================================
    //
    // Loop over low pt jet collection. Only for
    // MET Type-1 re-correction.
    //
    //==================================================
    for (int iJet=0; iJet < *nCorrT1METJet; iJet++){
      float recojet_pt_raw    = CorrT1METJet_rawPt[iJet];
      float recojet_eta       = CorrT1METJet_eta[iJet];
      float recojet_phi       = CorrT1METJet_phi[iJet];
      float recojet_area      = CorrT1METJet_area[iJet];
      float recojet_muonSubtrFactor = CorrT1METJet_muonSubtrFactor[iJet];

      jetCorrector->setJetPt(recojet_pt_raw);
      jetCorrector->setJetEta(recojet_eta);
      jetCorrector->setJetPhi(recojet_phi);
      jetCorrector->setJetA(recojet_area);
      jetCorrector->setRho(rho);
      float jec = jetCorrector->getCorrection();

      CorrT1METJet_pt[iJet] = jec * recojet_pt_raw;
      //========================================================================
      //
      // MET Type-1 correction
      //
      //========================================================================
      float recojet_cosphi = TMath::Cos(recojet_phi);
      float recojet_sinphi = TMath::Sin(recojet_phi);

      //
      // 1. Get the jet raw pt without muons included
      //
      float recojet_pt_noMuRaw   = recojet_pt_raw * (1.f - recojet_muonSubtrFactor);
      // This is not exactly the correct thing to do. We shouldn't be using the full
      // jet eta and phi but the eta and phi of muon-less jet p4. Not possible with
      // current NanoAODs.
      float recojet_eta_noMuRaw  = recojet_eta;
      float recojet_phi_noMuRaw  = recojet_phi;

      //
      // 2. Apply JEC on this muon-less jet raw p4
      // The JEC used is the same JEC we apply on the full raw p4
      float recojet_pt_noMuL1L2L3 = jec * recojet_pt_noMuRaw;
      float recojet_px_noMuL1L2L3 = recojet_pt_noMuRaw * TMath::Cos(recojet_phi_noMuRaw);
      float recojet_py_noMuL1L2L3 = recojet_pt_noMuRaw * TMath::Sin(recojet_phi_noMuRaw);
      CorrT1METJet_pt_noMuL1L2L3[iJet] = recojet_pt_noMuL1L2L3;

      //
      // 3.
      //
      float recojet_pt_noMuOnlyL1 = recojet_pt_noMuRaw;
      float recojet_px_noMuOnlyL1 = recojet_pt_noMuL1L2L3 * TMath::Cos(recojet_phi_noMuRaw);
      float recojet_py_noMuOnlyL1 = recojet_pt_noMuL1L2L3 * TMath::Sin(recojet_phi_noMuRaw);
      if (applyL1OnNoMuP4 && jetCorrectorL1){
        jetCorrectorL1->setJetPt(recojet_pt_raw);
        jetCorrectorL1->setJetEta(recojet_eta);
        jetCorrectorL1->setJetPhi(recojet_phi);
        jetCorrectorL1->setJetA(recojet_area);
        jetCorrectorL1->setRho(rho);
        float jecOnlyL1ForNoMuRaw = jetCorrectorL1->getCorrection();
        recojet_pt_noMuOnlyL1 *= jecOnlyL1ForNoMuRaw;
        recojet_px_noMuOnlyL1 = recojet_pt_noMuL1L2L3 * TMath::Cos(recojet_phi_noMuRaw);
        recojet_py_noMuOnlyL1 = recojet_pt_noMuL1L2L3 * TMath::Sin(recojet_phi_noMuRaw);
      }

      //
      // 4. Calculate the Type-1 correction and add it to the Raw MET
      //
      float recojet_emEF = 0.;
      if (hasEmEF){
        recojet_emEF = CorrT1METJet_EmEF[iJet];
      }

      if ((recojet_pt_noMuL1L2L3 > 15.f) && (recojet_emEF < 0.9f)){
        met_px_updated -= (recojet_px_noMuL1L2L3 - recojet_px_noMuOnlyL1);
        met_py_updated -= (recojet_py_noMuL1L2L3 - recojet_py_noMuOnlyL1);
      }
    }

    std::vector<float> v_jet_updated_idx_sortByPt(*nJet);
    std::iota(v_jet_updated_idx_sortByPt.begin(), v_jet_updated_idx_sortByPt.end(), 0);
    std::sort(v_jet_updated_idx_sortByPt.begin(), v_jet_updated_idx_sortByPt.end(),
      [&v_jet_pt_updated](float i1, float i2) {
        return v_jet_pt_updated[i1] > v_jet_pt_updated[i2];
      }
    );

    for (uint i=0; i < v_jet_updated_idx_sortByPt.size(); i++){
      Jet_IdxSortByPt[i] = v_jet_updated_idx_sortByPt[i];
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

  if(jetCorrector) delete jetCorrector;
  if(jetCorrectorL1) delete jetCorrectorL1;
  if(h2_jetvetomap) delete h2_jetvetomap;

  outTree->Write();
  outputTreeFile.Close();
  inputTreeFile.Close();
  std::cout << "AddJEC()::Finished" << std::endl;
}

