pcangsd="/home/users/xiaodong/Software/pcangsd/pcangsd.py"

indir="/home/users/xiaodong/Documents/Project/Seal/angsd/data"
outdir="/home/users/xiaodong/Documents/Project/Seal/pcangsd/data"

python ${pcangsd} -beagle $indir/world.clean.maf005.beagle.gz -o $outdir/test1 -threads 10 -iter 1000
