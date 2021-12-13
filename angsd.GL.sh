#!/usr/bin/bash

ANGSD=/home/users/xiaodong/Software/angsd/angsd
NGSADMIX=/home/users/xiaodong/Software/NGSadmixture/NGSadmix

# files
folder=/home/users/xiaodong/Documents/Project/Seal/angsd/data/20210120_r2


site=/home/users/xiaodong/Documents/Project/Seal/region_v2/20210120_r2/mappability.autosome.Newint_meantwo2_sd.subtractREPEAT.subtractF.r2.final.site
chr=/home/users/xiaodong/Documents/Project/Seal/region_v2/20210120_r2/mappability.autosome.Newint_meantwo2_sd.subtractREPEAT.subtractF.r2.final.chr

bamlist=/home/users/xiaodong/Documents/Project/Seal/sample_list_v2/20210119/r2/world.rmLow3andWSpecies.rmPlate1.v2.bamlist

bgl=${folder}/world.extra2.r2

# command
${ANGSD} -GL 1 -doMajorMinor 1 -doMaf 1 -SNP_pval 2e-6 -minMapQ 30 -minQ 20 -minInd 10 -minMaf 0.05 -doGlf 2\
	 -out $bgl -nThreads 10  -bam $bamlist -rf ${chr} -sites ${site} > ${bgl}.log


