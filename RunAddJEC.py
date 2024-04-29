import os
import glob
import argparse
import subprocess
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gSystem.Load("obj/libAddJEC.so")

inputFileDir="./output/"
outputFileDir="./output_JEC/"
jecVersion="Summer23BPixPrompt23_RunD_V1_DATA"

inFileList = [f for f in glob.glob(inputFileDir+"*")]
# inFileList = inFileList[0:1]

nInFiles = len(inFileList)

for iFile, inputFilePath in enumerate(inFileList):
  inputFileName=inputFilePath.split("/")[-1]
  outputFileName=inputFileName.replace("NanoSkim","NanoSkimJEC")
  outputFilePath=outputFileDir+outputFileName

  print(f"Run AddJEC() for file {iFile} / {nInFiles}")

  ROOT.AddJEC(inputFilePath,outputFilePath,"Events",jecVersion)
  print("")
