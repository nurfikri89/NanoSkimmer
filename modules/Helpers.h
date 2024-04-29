#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"

using namespace ROOT::VecOps;
using RNode = ROOT::RDF::RNode;

// RNode AddJEC(RNode &df){
//   std::string inDir("CondFormats/JetMETObjects/data");
//   std::string jecVersion("Summer23BPixPrompt23_RunD_V1_DATA");

//   JetCorrectorParameters Puppi_ResJetPar(inDir+"/"+jecVersion+"_L2L3Residual_AK4PFPuppi.txt");
//   JetCorrectorParameters Puppi_L3JetPar(inDir+"/"+jecVersion+"_L3Absolute_AK4PFPuppi.txt");
//   JetCorrectorParameters Puppi_L2JetPar(inDir+"/"+jecVersion+"_L2Relative_AK4PFPuppi.txt");
//   JetCorrectorParameters Puppi_L1JetPar(inDir+"/"+jecVersion+"_L1FastJet_AK4PFPuppi.txt");

//   // Load the JetCorrectorParameter objects into a vector, IMPORTANT: THE ORDER MATTERS HERE !!!!
//   std::vector<JetCorrectorParameters> vPar;
//   vPar.push_back(Puppi_L1JetPar);
//   vPar.push_back(Puppi_L2JetPar);
//   vPar.push_back(Puppi_L3JetPar);
//   vPar.push_back(Puppi_ResJetPar);
//   // FactorizedJetCorrector jetCorrector(vParAK4PFPuppi);

//   auto getJEC = [&vParTemp = vPar](RVec<float>& pt_raw, RVec<float>& eta, RVec<float>& phi, RVec<float>& area, RVec<float>& rho){
//     // thread_local auto jetCorrector = jc;
//     thread_local FactorizedJetCorrector jetCorrector(vParTemp);
//     RVec<float> pt_corr(pt_raw);
//     for (uint i=0; i<pt_raw.size();i++){
//       jetCorrector.setJetPt(pt_raw[i]);
//       jetCorrector.setJetEta(eta[i]);
//       jetCorrector.setJetPhi(phi[i]);
//       jetCorrector.setJetA(area[i]);
//       jetCorrector.setRho(rho[i]);
//       pt_corr[i] = pt_raw[i] * jetCorrector.getCorrection();
//     }
//     return pt_corr;
//   };
//   // define columns for hlt and l1 bitsets
//   df = df.Define("Jet_rho" ,    "RVec<float>(nJet,Rho_fixedGridRhoFastjetAll)");
//   df = df.Define("Jet_pt_raw" , "(1.f - Jet_rawFactor) * Jet_pt");

//   return df.Define("Jet_pt_nom", getJEC, {"Jet_pt_raw", "Jet_eta", "Jet_phi" , "Jet_area", "Jet_rho"});
// }


#include "nlohmann/json.hpp"
using json = nlohmann::json;
json golden_json;

void init_json(std::string jsonFile) {
  std::cout << "Initializing DC json file:" << jsonFile << std::endl;
  std::ifstream f(jsonFile);
  golden_json = json::parse(f);
}

bool isGoodLumi(int run, int lumi){
  for (auto& lumiRange : golden_json[std::to_string(run)]) {
    if (lumi >= lumiRange[0] && lumi <= lumiRange[1]) {
      return true;
    }
  }
  return false;
}