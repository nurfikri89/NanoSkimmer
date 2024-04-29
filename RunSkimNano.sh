echo ""
echo "=====START"
date
echo ""

SAMPLES=(
JetMET0_Run2024C-PromptReco-v1_NANOAOD
JetMET1_Run2024C-PromptReco-v1_NANOAOD
)

for SAMPLE in ${SAMPLES[@]};
do
  python3 -u SkimNano.py --sample $SAMPLE --cpus 4 --nfiles -1 2>&1 | tee logs/LOG_SkimNano_${SAMPLE}.txt
done

echo ""
echo ""
date
echo "=====END"
