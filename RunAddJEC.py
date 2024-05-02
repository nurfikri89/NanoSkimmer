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

# outputFileDir="./output_JEC_Winter23Prompt23/"
# jecVersion="Winter23Prompt23_RunA_V1_DATA"

inputFileDir="./nanoPrecision_input/"
outputFileDir="./nanoPrecision_output_Winter23_V2_AddJEC/"
jecVersion="Winter23Prompt23_V2_MC"

inFileList = [f for f in glob.glob(inputFileDir+"*.root")]
inFileList = inFileList[0:4]

nInFiles = len(inFileList)

for iFile, inputFilePath in enumerate(inFileList):
  inputFileName=inputFilePath.split("/")[-1]
  outputFileName=inputFileName.replace("NanoSkim","NanoSkimJEC")
  outputFilePath=outputFileDir+outputFileName

  print(f"Run AddJEC() for file {iFile} / {nInFiles}")

  ROOT.AddJEC(inputFilePath,outputFilePath,"Events",jecVersion)
  print("")

timerMain.Stop()
timerMain.Print()