import os
import glob
import argparse
import subprocess
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gSystem.Load("obj/libAddJEC.so")

timerMain = ROOT.TStopwatch()
timerMain.Start()

# inputFileDir="./output/"
# outputFileDir="./output_JEC_Summer23BPixPrompt23/"
# jecVersion="Summer23BPixPrompt23_RunD_V1_DATA"

inputFileDir="./output/"
outputFileDir="./output_JEC_Prompt24_V2M/"
jecVersion="Prompt24_V2M"

# inputFileDir="./nanoPrecision_input/"
# outputFileDir="./nanoPrecision_output_Winter23_V2_AddJEC/"
# jecVersion="Winter23Prompt23_V2_MC"

inFileList = [f for f in glob.glob(inputFileDir+"*.root")]
# inFileList = inFileList[0:1]
# inFileList = ["./output/NanoSkim_JetMET0_Run2024C-PromptReco-v1_NANOAOD_0201fc66-9eaa-445a-9a44-c1ced6cf5453.root"]

nInFiles = len(inFileList)

for iFile, inputFilePath in enumerate(inFileList):
  inputFileName=inputFilePath.split("/")[-1]
  isMC = "Summer" in inputFileName or "Winter" in inputFileName
  if isMC:
    print("This file is MC")
  else:
    print("This file is Data")
  outputFileName = inputFileName.replace("NanoSkim","NanoSkimJEC")
  outputFilePath=outputFileDir+outputFileName

  print(f"Run AddJEC() for file {iFile} / {nInFiles}")

  ROOT.AddJEC(inputFilePath,outputFilePath,"Events",jecVersion,isMC)
  print("")

timerMain.Stop()
timerMain.Print()