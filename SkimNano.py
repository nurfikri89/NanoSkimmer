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
parser.add_argument('-m', '--isMC',     action='store_true', default=False)
args = parser.parse_args()

sampleName = args.sample
nfiles = args.nfiles
isMC = args.isMC

# goldenJSONPath  = "data/lumi/Collisions24_13p6TeV_378981_380074_DCSOnly_TkPx.json"
# goldenJSONPath  = "data/lumi/Cert_Collisions2024_378981_379470_Golden.json"
# goldenJSONPath  = "data/lumi/Cert_Collisions2024_378981_379866_Golden.json"
# goldenJSONPath  = "data/lumi/Cert_Collisions2024_378981_380115_Golden.json"
goldenJSONPath  = "data/lumi/Cert_Collisions2024_378981_380470_Golden.json"

inFiles=[]
with open(f"./samples/{sampleName}.txt", 'r') as txtfile:
  inFiles = [line.rstrip() for line in txtfile]

if nfiles > 0:
  print(f"Only process {nfiles} files")
  inFiles = inFiles[0:nfiles]
else:
  print("Process all files")

PathPrefix=""
inFiles = [f for f in inFiles if not(f.startswith("#"))]
inFilesFinal = [PathPrefix+f for f in inFiles]
nFiles = len(inFilesFinal)

print(f"Processing sample {sampleName}")
print(f"isMC = {isMC}")

##########################################
#
#
#
##########################################
ROOT.ROOT.EnableImplicitMT(args.cpus)

if not(isMC):
  print(f"Setup DC JSON: {goldenJSONPath}")
  ROOT.init_json(goldenJSONPath)

def SkimNanoFile(nFiles, iFile, inFilePath):
  inFileName = inFilePath.split('/')[-1]
  inFileName = inFileName.replace(".root","")

  print("")
  print("==============================================")
  print(f"Skimming {iFile}/{nFiles} : {inFilePath}")

  df = ROOT.ROOT.RDataFrame("Events",inFilePath)
  # ROOT.RDF.Experimental.AddProgressBar(df) #Only available in ROOT version >= 6.30
  df_count_initial = df.Count()

  if not(isMC):
    print("Filter on Golden JSON")
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
    # isHLTNotJetFlag = ("HLT_" in bStr or "DST_" in bStr) and not("PFJet" in bStr)
    isHLTNotJetFlag = ("HLT_" in bStr or "DST_" in bStr) and not("HLT_PFJet" in bStr)
    isL1Flag = "L1_" in bStr
    isLowPtElec = "nLowPtElectron" in bStr or "LowPtElectron_" in bStr
    isSV = "nSV" in bStr or "SV_" in bStr
    isSoftActivityJet = "SoftActivityJet" in bStr
    isSubJet = "nSubJet" in bStr or "SubJet_" in bStr
    isFatJet = "nFatJet" in bStr or "FatJet_" in bStr
    isTrigObj = "nTrigObj" in bStr or "TrigObj_" in bStr
    isIsoTrack = "nIsoTrack" in bStr or "IsoTrack_" in bStr
    isPPS = "nProton" in bStr or "Proton_" in bStr
    isTauProd = "nTauProd" in bStr or "TauProd_" in bStr
    isTau = "nTau" in bStr or "Tau_" in bStr
    isBoostedTau = "nboostedTau" in bStr or "boostedTau_" in bStr
    isHTXS = "HTXS" in bStr
    isGenPart = "nGenPart" in bStr or "GenPart_" in bStr
    isGenJetAK8 = "nGenJetAK8" in bStr or "GenJetAK8_" in bStr
    isSubGenJet = "nSubGenJetAK8" in bStr or "SubGenJetAK8_" in bStr
    isFsrPhoton = "nFsrPhoton" in bStr or "FsrPhoton_" in bStr
    isGenDressedLepton = "nGenDressedLepton" in bStr or "GenDressedLepton_" in bStr
    isGenIsolatedPhoton = "nGenIsolatedPhoton" in bStr or "GenIsolatedPhoton_" in bStr
    isJetCHS = "nJetCHS" in bStr or "JetCHS_" in bStr
    isJetMisc = "Jet_puId" in bStr or "Jet_qgl" in bStr or "Jet_svIdx" in bStr
    isFatJetForJEC =  "nFatJetForJEC" in bStr or "FatJetForJEC_" in bStr
    isGenJetAK8ForJEC = "nGenJetAK8ForJEC" in bStr or "GenJetAK8ForJEC_" in bStr
    isJetCalo = "nJetCalo" in bStr or "JetCalo_" in bStr
    if not(isHLTNotJetFlag or isL1Flag or isLowPtElec or isTrigObj or \
      isSV or isSoftActivityJet or isSubJet or isFatJet or isTau or isBoostedTau or \
      isFsrPhoton or isGenDressedLepton or isGenIsolatedPhoton or \
      isGenPart or isGenJetAK8 or \
      isSubGenJet or isIsoTrack or isJetCHS or isJetMisc or \
      isFatJetForJEC or isGenJetAK8ForJEC or \
      isJetCalo or isPPS or isTauProd or isHTXS):
      # print(f"{bStr}")
      branchesToSaveFinal.append(bStr)
  rdf_opts = ROOT.RDF.RSnapshotOptions()
  rdf_opts.fLazy = True
  rdf_opts.fCompressionLevel = 9
  rdf_opts.fCompressionAlgorithm = 2
  df = df.Snapshot("Events", outFilePath, branchesToSaveFinal, rdf_opts)
  # print("After snapshot")

  ncount_initial=df_count_initial.GetValue()
  ncount_final=df_count_final.GetValue()

  print(f"{ncount_initial=}")
  print(f"{ncount_final=}")

  if ncount_final == 0:
    print(f"No events at the end. Deleting output file: {outFilePath}")
    os.remove(outFilePath)

  del df

print("\n")
print(f"Skimming nFiles = {nFiles}")
for iFile, inFile in enumerate(inFilesFinal):
  SkimNanoFile(nFiles, iFile, inFile)
