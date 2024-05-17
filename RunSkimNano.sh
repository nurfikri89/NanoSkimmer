echo ""
echo "=====START"
date
echo ""

SAMPLES=(
)
for SAMPLE in ${SAMPLES[@]};
do
  python3 -u SkimNano.py --sample $SAMPLE --cpus 4 --nfiles 2 2>&1 | tee logs/LOG_SkimNano_${SAMPLE}.txt
done


SAMPLES=(
QCD_Flat2022_Run3Winter24JMENano_GTv10
QCD_Flat2022_Run3Winter24JMENano_GTv9
)

for SAMPLE in ${SAMPLES[@]};
do
  python3 -u SkimNano.py --isMC --sample $SAMPLE --cpus 4 --nfiles -1 2>&1 | tee logs/LOG_SkimNano_${SAMPLE}.txt
done


echo ""
echo ""
date
echo "=====END"
