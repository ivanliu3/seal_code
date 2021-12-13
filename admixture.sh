#!/usr/bin/bash

K=$1 # K given by argument
ANGSD=/home/users/xiaodong/Software/angsd/angsd
NGSADMIX=/home/users/xiaodong/Software/NGSadmixture/NGSadmix

# files
angsdfolder=/home/users/xiaodong/Documents/Project/Seal/angsd/data/20210120_r2
admixfolder=/home/users/xiaodong/Documents/Project/Seal/admixture/data/2021_r2/world

# command

for i in {1..200}
do
    echo $i >> $admixfolder/r2.extra2base.world.${K}.log
    $NGSADMIX -likes ${angsdfolder}/world.extra2.r2.beagle.gz  -K ${K} -P 4 -maxiter 5000 -minMaf 0.05  -o $admixfolder/r2.extra2base.world.K${K}.round.$i &>> $admixfolder/r2.extra2base.world.${K}.log
    echo >> $admixfolder/r2.extra2base.world.${K}.log
done






