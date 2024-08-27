import os
import glob
import argparse
import subprocess
import ROOT
ROOT.gROOT.SetBatch()
ROOT.gSystem.Load("obj/libHelpers.so")

parser = argparse.ArgumentParser("")
parser.add_argument('-s', '--sample',   type=str,  default="")
parser.add_argument('-c', '--cpus',     type=int,  default=1)
parser.add_argument('-f', '--nfiles',   type=int,  default=-1)
parser.add_argument('-m', '--isMC',     action='store_true', default=False)
parser.add_argument('--doSkimTeVJet',   action='store_true', default=False)
parser.add_argument('--doSkimPhoton',   action='store_true', default=False)
parser.add_argument('--doSkimDiElec',   action='store_true', default=False)

args = parser.parse_args()

downloadInput = False

sampleName = args.sample
nfiles = args.nfiles
isMC = args.isMC

doSkimTeVJet = args.doSkimTeVJet
doSkimPhoton = args.doSkimPhoton
doSkimDiElec = args.doSkimDiElec

goldenJSONPath2024E = "data/lumi/Collisions24_13p6TeV_378981_381309_DCSOnly_TkPx.json"
goldenJSONPath2024  = "data/lumi/Cert_Collisions2024_378981_380470_Golden.json"
goldenJSONPath2023  = "data/lumi/Cert_Collisions2023_366442_370790_Golden.json"

# lxplus
# outDir = f"/eos/cms/store/group/phys_jetmet/nbinnorj/NanoSkimJEC/"
# TMPDIR=os.getenv("TMPDIR")

#Hefaistos
TMPDIR="/tmp/"

applyTeVJetSel = False
applyPhoton200Sel = False
applyDiElectron30 = False

if doSkimTeVJet:
  # outDir = f"./output_JetTeVSkim/"
  # outDir = f"./output_JetTeVSkim_PFNano/"
  outDir = f"./output_JetTeVSkim_PFNano_0p2/"
  applyTeVJetSel = True
if doSkimPhoton:
  outDir = f"./output_Photon200Skim/"
  applyPhoton200Sel = True

if doSkimDiElec:
  outDir = f"./output_DiElectron30JetPt180Skim/"
  applyDiElectron30 = True

if not(applyTeVJetSel) and not(applyPhoton200Sel) and not(applyDiElectron30):
  raise Exception(f"No selection is applied. Please check! {applyTeVJetSel=}, {applyPhoton200Sel=}, {applyDiElectron30=}")


outDirFinal = f"{outDir}/{sampleName}"
if not os.path.exists(outDirFinal):
  os.makedirs(outDirFinal)
##########################################
#
#
#
##########################################
# inFilesDir="./samples"
inFilesDir="./samplesPFNano/"

inFiles=[]
with open(f"{inFilesDir}/{sampleName}.txt", 'r') as txtfile:
  inFiles = [line.rstrip() for line in txtfile]

if nfiles > 0:
  print(f"Only process {nfiles} files")
  inFiles = inFiles[0:nfiles]
else:
  print("Process all files")

#
#
#
inFilesFinal = [f for f in inFiles if not(f.startswith("#"))]
nFiles = len(inFilesFinal)

##########################################
#
#
#
##########################################
print(f"Processing sample {sampleName}")
print(f"isMC = {isMC}")
ROOT.ROOT.EnableImplicitMT(args.cpus)

goldenJSONPath = None
if "Run2024" in sampleName:
  goldenJSONPath = goldenJSONPath2024
if "Run2024E" in sampleName:
  goldenJSONPath = goldenJSONPath2024E
if "Run2023" in sampleName:
  goldenJSONPath = goldenJSONPath2023

if not(isMC):
  print(f"Setup DC JSON: {goldenJSONPath}")
  ROOT.init_json(goldenJSONPath)

def SkimNanoFile(nFiles, iFile, inFilePath):
  print("")
  print("==============================================")
  print(f"Skimming {iFile}/{nFiles} : {inFilePath}")

  inFileName = inFilePath.split('/')[-1]
  inFileName = inFileName.replace(".root","")

  PathPrefix="root://xrootd-cms.infn.it/"

  inFileInTMPDIR = ""
  inFilePathFinal = None

  if inFilePath.startswith("/store/"):
    inFilePathTmp = PathPrefix+inFilePath
    if downloadInput:
      inFilePathFinal = f"{TMPDIR}/{inFileName}.root"
      print(f"Copying input file {inFilePathTmp}")
      print(f"to tmp dir: {inFilePathFinal}")
      command = ["xrdcp", "-f", inFilePathTmp, inFilePathFinal]
      subprocess.run(command)
      inFileInTMPDIR = inFilePathFinal
    else:
      inFilePathFinal = inFilePathTmp
  else:
    inFilePathFinal = inFilePath

  #
  #
  #
  tempInfile = ROOT.TFile(inFilePathFinal,"OPEN")
  tempInTree = tempInfile.Get("Events")
  tempNEntries = tempInTree.GetEntries()

  if tempNEntries == 0:
    print(f"Events tree in {inFilePathFinal} has 0 entries. Skipping this file")
    return

  df = ROOT.ROOT.RDataFrame("Events",inFilePathFinal)
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
  TriggerFlags = []
  if applyTeVJetSel:
    TriggerFlags += ["HLT_PFJet500"]
    df = df.Define("nJetPt1000","Sum(Jet_pt >= 1000.f)")
    df = df.Define("passTeVJetSel","nJetPt1000>=1")

  if applyPhoton200Sel:
    TriggerFlags += [
      # "HLT_Photon50EB_TightID_TightIso",
      # "HLT_Photon55EB_TightID_TightIso",
      # "HLT_Photon75EB_TightID_TightIso",
      # "HLT_Photon75EB_TightID_TightIso",
      # "HLT_Photon90EB_TightID_TightIso",
      "HLT_Photon200"
    ]
    df = df.Define("nTightPhoton200","Sum(Photon_pt >= 200.f && Photon_cutBased >= 3)")
    df = df.Define("passPhotonSel","nTightPhoton200>=1")

  if applyDiElectron30:
    TriggerFlags += ["HLT_Ele30_WPTight_Gsf","HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL"]
    df = df.Define("nTightElectron30","Sum(Electron_pt >= 30.f && Electron_cutBased >= 3)")
    df = df.Define("nJetPt180","Sum(Jet_pt >= 180.f)")
    df = df.Define("passDiElectronSel","nTightElectron30>=2 && nJetPt180>=1")

  if len(TriggerFlags) > 0:
    TriggerFlagsSel = "("+"||".join(TriggerFlags)+")"
  else:
    TriggerFlagsSel = "(true)"

  #
  #
  #
  SelectionStr = f"{TriggerFlagsSel} && {METFiltersSel} && nJet>=1"
  if applyTeVJetSel:
    SelectionStr += " && passTeVJetSel"
  if applyPhoton200Sel:
    SelectionStr += " && passPhotonSel"
  if applyDiElectron30:
    SelectionStr += " && passDiElectronSel"

  print(f"Applying selection : {SelectionStr}")
  df = df.Filter(SelectionStr)

  # df = ROOT.AddJEC(ROOT.RDF.AsRNode(df))
  df_count_final = df.Count()

  #
  #
  #
  outFileName=f"NanoSkim_{sampleName}_{inFileName}.root"
  outFilePathTemp=f"{TMPDIR}/{outFileName}"

  outFilePathFinal=""
  if "/eos/user" in outDirFinal:  outFilePathFinal+="root://eosuser.cern.ch/"
  elif "/eos/cms" in outDirFinal: outFilePathFinal+="root://eoscms.cern.ch/"
  outFilePathFinal += f"{outDirFinal}/{outFileName}"

  branchesToSave = df.GetColumnNames()
  branchesToSaveFinal = []
  for b in branchesToSave:
    bStr = str(b)
    isHLTFlag = False
    if applyTeVJetSel:
      isHLTFlag = ("HLT_" in bStr or "DST_" in bStr) and not("HLT_PFJet" in bStr)
    if applyDiElectron30:
      isHLTFlag = ("HLT_" in bStr or "DST_" in bStr) and not("HLT_Ele" in bStr)
    if applyPhoton200Sel:
      isHLTFlag = ("HLT_" in bStr or "DST_" in bStr) and not("HLT_Photon" in bStr)
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
    if not(isHLTFlag or isL1Flag or isLowPtElec or isTrigObj or \
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
  df = df.Snapshot("Events", outFilePathTemp, branchesToSaveFinal, rdf_opts)
  # print("After snapshot")

  ncount_initial=df_count_initial.GetValue()
  ncount_final=df_count_final.GetValue()

  print(f"{ncount_initial=}")
  print(f"{ncount_final=}")

  if ncount_final == 0:
    print(f"No events at the end. Deleting output file: {outFilePathTemp}")
    os.remove(outFilePathTemp)
  else:
    print(f"Copying output file to final path: {outFilePathFinal}")
    command = ["xrdcp", "-f", outFilePathTemp, outFilePathFinal]
    subprocess.run(command)
    print(f"Deleting output file in TMPDIR: {outFilePathTemp}")
    os.remove(outFilePathTemp)

  if inFileInTMPDIR != "":
    print(f"Deleting input file in TMPDIR: {inFileInTMPDIR}")
    os.remove(inFileInTMPDIR)

  del df

print("\n")
print(f"Skimming nFiles = {nFiles}")
for iFile, inFile in enumerate(inFilesFinal):
  SkimNanoFile(nFiles, iFile, inFile)
