import os
import glob
import argparse
import subprocess
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gSystem.Load("obj/libAddJEC_ForNanoValidation.so")

timerMain = ROOT.TStopwatch()
timerMain.Start()

inputFileDir="./nanoPrecision_input/"
outputFileDir="./nanoPrecision_output_Winter23_V2_AddJEC_ForNanoValidation/"
jecVersion="Winter23Prompt23_V2_MC"

inFileList = [f for f in glob.glob(inputFileDir+"*.root")]
# inFileList = inFileList[0:4]

nInFiles = len(inFileList)

for iFile, inputFilePath in enumerate(inFileList):
  inputFileName=inputFilePath.split("/")[-1]
  outputFileName=inputFileName.replace("NanoSkim","NanoSkimJEC")
  outputFileName=inputFileName.replace("NanoAOD","NanoAODJEC")
  outputFilePath=outputFileDir+outputFileName

  print(f"Run AddJEC_ForNanoValidation() for file {iFile} / {nInFiles}")

  ROOT.AddJEC_ForNanoValidation(inputFilePath,outputFilePath,"Events",jecVersion)
  print("")

timerMain.Stop()
timerMain.Print()