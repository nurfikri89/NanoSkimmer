import os
import tarfile
import urllib.request

import urllib.request

baseURL="https://raw.githubusercontent.com/cms-jet/JECDatabase/master/tarballs/"
jesFilesList=[
  "Summer22_22Sep2023_V2_MC.tar.gz",
  "Summer22EE_22Sep2023_V2_MC.tar.gz",
  "Summer23Prompt23_V1_MC.tar.gz",
  "Summer23BPixPrompt23_V1_MC.tar.gz",
  "Summer23BPixPrompt23_RunD_V1_DATA.tar.gz"
]

for jesFile in jesFilesList:
  print(f"Downloading {baseURL}/{jesFile}")
  urllib.request.urlretrieve(f"{baseURL}/{jesFile}",jesFile)
  print(f"Untar {jesFile}")
  jesArchive = tarfile.open(jesFile, "r:gz")
  jesArchive.extractall("./")
  print(f"Deleting {jesFile}")
  os.remove(jesFile)
