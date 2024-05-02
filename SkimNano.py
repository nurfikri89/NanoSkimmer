import os
import glob
import argparse
import ROOT
ROOT.gROOT.SetBatch()
ROOT.gSystem.Load("obj/libHelpers.so")

XROOTD="root://xrootd-cms.infn.it/"

parser = argparse.ArgumentParser("")
parser.add_argument('-s', '--sample',   type=str,  default="")
parser.add_argument('-c', '--cpus',     type=int,  default=1)
parser.add_argument('-f', '--nfiles',   type=int,  default=-1)
args = parser.parse_args()

sampleName = args.sample
nfiles = args.nfiles

# goldenJSONPath  = "data/lumi/Collisions24_13p6TeV_378981_380074_DCSOnly_TkPx.json"
goldenJSONPath  = "data/lumi/Cert_Collisions2024_378981_379470_Golden.json"

inFiles=[]
with open(f"./samples/{sampleName}.txt", 'r') as txtfile:
  inFiles = [line.rstrip() for line in txtfile]

if nfiles > 0:
  print(f"Only process {nfiles} files")
  inFiles = inFiles[0:nfiles]
else:
  print("Process all files")

PathPrefix=""
inFilesFinal = [PathPrefix+f for f in inFiles]
nFiles = len(inFilesFinal)

print(f"Processing sample {sampleName}")

##########################################
#
#
#
##########################################
ROOT.ROOT.EnableImplicitMT(args.cpus)

def SkimNanoFile(nFiles, iFile, inFilePath):
  inFileName = inFilePath.split('/')[-1]
  inFileName = inFileName.replace(".root","")

  print("")
  print("==============================================")
  print(f"Skimming {iFile}/{nFiles} : {inFilePath}")

  df = ROOT.ROOT.RDataFrame("Events",inFilePath)
  # ROOT.RDF.Experimental.AddProgressBar(df) #Only available in ROOT version >= 6.30
  df_count_initial = df.Count()

  #
  #
  #
  ROOT.init_json(goldenJSONPath)
  df = df.Define("passGoldenJSON", "isGoodLumi(run, luminosityBlock)")
  df = df.Filter("passGoldenJSON")

  #
  #
  #
  METFilters = [
    "Flag_goodVertices",
    "Flag_globalSuperTightHalo2016Filter",
    "Flag_EcalDeadCellTriggerPrimitiveFilter",
    "Flag_BadPFMuonFilter",
    "Flag_BadPFMuonDzFilter",
    "Flag_hfNoisyHitsFilter",
    "Flag_eeBadScFilter",
    "Flag_ecalBadCalibFilter",
  ]
  METFiltersSel = "("+"&&".join(METFilters)+")"

  #
  #
  #
  TriggerFlags = [
    "HLT_PFJet500",
  ]
  TriggerFlagsSel = "("+"||".join(TriggerFlags)+")"


  #
  #
  #
  SelectionStr = f"{TriggerFlagsSel} && {METFiltersSel} && nJet>=1"
  df = df.Filter(SelectionStr)

  # df = ROOT.AddJEC(ROOT.RDF.AsRNode(df))
  df_count_final = df.Count()

  outFilePath=f"./output/NanoSkim_{sampleName}_{inFileName}.root"
  branchesToSave = df.GetColumnNames()

  branchesToSaveFinal = []
  for b in branchesToSave:
    bStr = str(b)
    isHLTNotJetFlag = ("HLT_" in bStr or "DST_" in bStr) and not("PFJet" in bStr)
    isL1Flag =  "L1_" in bStr
    isLowPtElec =  "LowPtElectron_" in bStr
    isSV =  "nSV" in bStr or "SV_" in bStr
    isSubJet =  "nSubJet" in bStr or "SubJet_" in bStr
    isFatJet =  "nFatJet" in bStr or "FatJet_" in bStr
    isPPS =  "nProton" in bStr or "Proton_" in bStr
    isTauProd =  "nTauProd" in bStr or "TauProd_" in bStr
    if not(isHLTNotJetFlag or isL1Flag or isLowPtElec or \
      isSV or isSubJet or isFatJet or isPPS or isTauProd):
      branchesToSaveFinal.append(bStr)
  rdf_opts = ROOT.RDF.RSnapshotOptions()
  rdf_opts.fLazy = True
  rdf_opts.fCompressionLevel = 9
  rdf_opts.fCompressionAlgorithm = 2
  df = df.Snapshot("Events", outFilePath, branchesToSaveFinal, rdf_opts)
  # print("After snapshot")

  print(f"{df_count_initial.GetValue()=}")
  print(f"{df_count_final.GetValue()=}")

  del df



print("\n")
print(f"Skimming nFiles = {nFiles}")
for iFile, inFile in enumerate(inFilesFinal):
  SkimNanoFile(nFiles, iFile, inFile)
