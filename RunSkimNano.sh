echo ""
echo "=====START"
date
echo ""

SAMPLES=(
# JetMET0_Run2024C-PromptReco-v1_NANOAOD
# JetMET0_Run2024D-PromptReco-v1_NANOAOD
# JetMET0_Run2023D-22Sep2023_v1-v1_NANOAOD
# JetMET0_Run2023D-22Sep2023_v2-v1_NANOAOD
# JetMET0_Run2023D-19Dec2023-v1_NANOAOD
)
NFILES=-1

for SAMPLE in ${SAMPLES[@]};
do
  python3 -u SkimNano.py --doSkimTeVJet --sample $SAMPLE --cpus 4 --nfiles ${NFILES} 2>&1 | tee logs/LOG_SkimNano_${SAMPLE}.txt
done

SAMPLES=(
EGamma0_Run2024C
EGamma0_Run2024D
EGamma0_Run2023D-22Sep2023_v1-v1_NANOAOD
EGamma0_Run2023D-22Sep2023_v2-v1_NANOAOD
)
NFILES=-1

for SAMPLE in ${SAMPLES[@]};
do
  python3 -u SkimNano.py --doSkimPhoton --sample $SAMPLE --cpus 4 --nfiles ${NFILES} 2>&1 | tee logs/LOG_SkimNano_${SAMPLE}.txt
done

SAMPLES=(
# EGamma0_Run2024C
# EGamma0_Run2024D
# EGamma0_Run2023D-22Sep2023_v1-v1_NANOAOD
# EGamma0_Run2023D-22Sep2023_v2-v1_NANOAOD
)
NFILES=-1

for SAMPLE in ${SAMPLES[@]};
do
  python3 -u SkimNano.py --doSkimDiElec --sample $SAMPLE --cpus 4 --nfiles ${NFILES} 2>&1 | tee logs/LOG_SkimNano_${SAMPLE}.txt
done

echo ""
echo ""
date
echo "=====END"
