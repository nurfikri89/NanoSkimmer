import os
import ROOT
ROOT.gROOT.SetBatch(True)

cwd = os.getcwd()

ROOT.gROOT.ProcessLine(".O 2") # Set optimization level to 2
ROOT.gROOT.ProcessLine(".O") # Show optimization level
ROOT.gSystem.SetBuildDir(f"{cwd}/obj/",True);

ROOT.gSystem.AddIncludePath(f" -I{cwd}/")
ROOT.gSystem.AddIncludePath(f" -I{cwd}/json/include/")

# ROOT.gROOT.ProcessLine(f".L {cwd}/CondFormats/JetMETObjects/src/Utilities.cc+O");
# ROOT.gROOT.ProcessLine(f".L {cwd}/CondFormats/JetMETObjects/src/JetCorrectorParameters.cc+O");
# ROOT.gROOT.ProcessLine(f".L {cwd}/CondFormats/JetMETObjects/src/SimpleJetCorrector.cc+O");
# ROOT.gROOT.ProcessLine(f".L {cwd}/CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc+O");
# ROOT.gROOT.ProcessLine(f".L {cwd}/src/Helpers.h+O")
# ROOT.gROOT.ProcessLine(f".L {cwd}/src/AddJEC.C+O")

ROOT.gSystem.CompileMacro(f"{cwd}/CondFormats/JetMETObjects/src/Utilities.cc","kO+","libUtilities");
ROOT.gSystem.CompileMacro(f"{cwd}/CondFormats/JetMETObjects/src/JetCorrectorParameters.cc","kO+","libJetCorrectorParameters");
ROOT.gSystem.CompileMacro(f"{cwd}/CondFormats/JetMETObjects/src/SimpleJetCorrector.cc","kO+","libSimpleJetCorrector");
ROOT.gSystem.CompileMacro(f"{cwd}/CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc","kO+","libFactorizedJetCorrector");
ROOT.gSystem.CompileMacro(f"{cwd}/modules/Helpers.h","kO+","libHelpers");
ROOT.gSystem.CompileMacro(f"{cwd}/modules/AddJEC.C","kO+","libAddJEC");
