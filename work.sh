# reads mapping
bwa mem -t 4 -k 32 -M -R "@RG\tID:ReadsGroup\tLB:Library\tSM:SampleID" reference.fa read1.fq.gz read2.fq.gz | samtools view -b -t reference.fa.fai - > raw.bam
samtools sort -m 4000000000 raw.bam sorted.bam
samtools rmdup sorted.bam rmdup.bam

# snp calling
java -jar -Xmx4g GenomeAnalysisTK.jar -T HaplotypeCaller --variant_index_type LINEAR --variant_index_parameter 128000 -R reference.fa -I rmdup.bam --emitRefConfidence GVCF -o raw.g.vcf
java -jar -Xmx7g GenomeAnalysisTK.jar -T GenotypeGVCFs -R reference.fa -V gvcf.list -o raw.vcf

# snp filter
java -jar -Xmx4g GenomeAnalysisTK.jar -T SelectVariants -R reference.fa -V raw.vcf -selectType SNP -o snp.vcf
java -jar -Xmx4g GenomeAnalysisTK.jar -T VariantFiltration -R reference.fa -V snp.vcf --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filterName "my_snp_filter" -o snp.filtered.vcf

# admixture
for i in {2..15}
do
	admixture snp.bed $i --cv -j4 -C 0.000001
done

# pca
flashpca --bfile snp --numthreads 8 --suffix _suffix

# tree
treebest nj -b 100 aligned.snp.fa > nj_tree.out

# snp imputation and phase
java -Xss5m -Xmx10g -jar beagle.jar gtgl=filtered.vcf out=imputed.vcf nthreads=16
java -Xss5m -Xmx10g -jar beagle.jar gt=imputed.vcf out=phased.vcf nthreads=16

# ld block
plink --noweb --bfile snp --blocks no-pheno-req

# identity score (IS) was calculated by an in-house python scripths
python cal_ISmatrix_2.0.py snp.vcf snp.vcf.index output

# pair wise IBS similarity
plink --noweb --bfile snp --genome --out output

# nucleotide diversity
vcftools --vcf subpop.vcf --window-pi 20000 --window-pi-step 10000 --out output

# Fst
vcftools --vcf snp.vcf --weir-fst-pop subpop1.lst --weir-fst-pop subpop2.lst --fst-window-size 20000 --fst-window-step 10000 --out output

# TajimaD 
vcftools --vcf snp.vcf --TajimaD --out output

# ROH
plink --bfile snp --homozyg-density 50 --homozyg-window-het 0 --homozyg-snp 50 --homozyg-window-snp 50 --homozyg-window-missing 2 --homozyg-kb 1 --homozyg-gap 5000

# ld
java -Xmx12g -jar Haploview.jar -memory 2000 -maxdistance 500 -n -log log.txt -pedfile pop.ped -info pop.info -out output -dprime -minMAF  0.05

# treemix
treemix -i input.gz -o output -k 100 -m 0 -root CA -global

# introgression(fd)
python egglib_sliding_windows.py -i input -o output.csv \
-w 20000 -m 60 -s 10000 --report 100 -a pi,dxy,ABBABABA,popS,S \
-p "P1[NW02,PW05,PW06];P2[NP01,NP03,PL01,PL02,PL03,PL04,PL05,PL06];P3[GHT-1,GHT-2,GHT-3,GHT-4,PW08];O[DL01,DL02,DL03,DL04,DL05,DL06,DL07,DL08,DL09,DL10,DL11,DL12,DB13,DB14,DB15];P3a[GHT-1,GHT-2];P3b[GHT-3,GHT-4,PW08]" \
--minimumExploitableData 0.5

# ihs
selscan --ihs --vcf input.vcf --map input.map --out outprefix --threads 4 --maf 0.05 --skip-low-freq
norm --ihs --files outprefix* --bp-win --winsize 20000

# nsl
selscan --nsl --vcf input.vcf --map input.map --out outprefix --threads 4 --maf 0.05 --skip-low-freq
norm --ihs --files outprefix* --bp-win --winsize 20000

# xpehh
selscan --xpehh --hap input.hap --ref ref.hap --map input.map --maf 0.05 --skip-low-freq --out output --threads 4
norm --xpehh --files outprefix* --bp-win --winsize 20000

# xpclr
XPCLR -xpclr pop1.geno pop2.geno snp.map outprefix -w1 0.0005 200 200 1 -p1 0.9

