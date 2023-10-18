while getopts :r:p:g:o: opt
do
case ${opt} in
r) ROBUST=${OPTARG};;
p) PERMISSIVE=${OPTARG};;
g) REFER_GFF=${OPTARG};;
o) OUTDIR=${OPTARG};;
*) echo "Usage is : $0 -r <robust peak> -p <permissive peak> -g <reference gff> -o <output directory> ">&2
    exit 1;;
esac
done

# step1 GFF Araport 11 
mkdir $OUTDIR/01.refer_stat
cd $OUTDIR/01.refer_stat
>gff_annotation.txt
for feature in gene mRNA five_prime_UTR;
do
    echo $feature `cat ${REFER_GFF} | awk '$3=="'$feature'"' | wc -l` >> gff_annotation.txt
done 
cat ${REFER_GFF} | awk '$3=="gene"' > Araport11_GFF3_gene
cat ${REFER_GFF} | awk '$3=="mRNA"' > Araport11_GFF3_mRNA
cat ${REFER_GFF} | awk '$3=="five_prime_UTR"' > Araport11_GFF3_five_prime_UTR

# step2 BED to GFF
mkdir $OUTDIR/02.bed2gff
cd $OUTDIR/02.bed2gff
cp $ROBUST .
gunzip tc.spi_merged.ctssMaxCounts11_ctssMaxTpm1.bed.gz
gt bed_to_gff3 -force -o tc.spi_merged.ctssMaxCounts11_ctssMaxTpm1.gff3 tc.spi_merged.ctssMaxCounts11_ctssMaxTpm1.bed

# step3 CAGE和GFF交集
mkdir $OUTDIR/03.cage_gff
cd $OUTDIR/03.cage_gff
## setp 3.1 Robust peak
bedtools intersect -a $OUTDIR/02.bed2gff/tc.spi_merged.ctssMaxCounts11_ctssMaxTpm1.gff3 -b ${REFER_GFF} -wb > bbb
## setp3.2 交集统计
>cage_gff_intersect.txt
for feature in gene mRNA five_prime_UTR;
do 
   echo $feature `cat bbb | awk '$12=="'$feature'"' | awk '{print $18}' | awk -F';' '{print $1}' | sed 's/ID=//' | sort | uniq | wc -l` >> cage_gff_intersect.txt
done

# step4 peak 统计
mkdir $OUTDIR/04.robust
cd $OUTDIR/04.robust
## step4.1 peak in utr5 gene numbers 
cat $OUTDIR/03.cage_gff/bbb | awk '$12=="five_prime_UTR"' | awk '{print $18}' | awk -F';' '{print $1}' | sed 's/ID=//' | sort | uniq | awk -F':' '{print $1}' | sort | uniq | wc -l 
## step4.2 coding gene or nonconding gene 
cat $OUTDIR/03.cage_gff/bbb | awk '$12=="five_prime_UTR"' | awk '{print $18}' | awk -F';' '{print $1}' | sed 's/ID=//' | sort | uniq | awk -F':' '{print $1}' | sort | uniq > robust_utr5_genelist
cat $OUTDIR/01.refer_stat/Araport11_GFF3_gene | awk '{split($9,a,";");split(a[1],b,"=");{print b[2]"\t"$0}}' > Araport11_GFF3_gene_addid

## step4.3 peaks 提取 BED_feature
cat $OUTDIR/02.bed2gff/tc.spi_merged.ctssMaxCounts11_ctssMaxTpm1.gff3 | awk '$3=="BED_feature"' > robust_BED_feature.gff3

## step4.4 BED_feature 和 5UTR区间作为p1@genename
bedtools intersect -a $OUTDIR/01.refer_stat/Araport11_GFF3_five_prime_UTR -b $OUTDIR/02.bed2gff/tc.spi_merged.ctssMaxCounts11_ctssMaxTpm1.bed -wb > p_format_in_fantom_prepare.csv

## setp4.5 非5' utr区间 p@locus
bedtools intersect -a $OUTDIR/01.refer_stat/Araport11_GFF3_mRNA -b $OUTDIR/02.bed2gff/tc.spi_merged.ctssMaxCounts11_ctssMaxTpm1.bed -wb > p_format_in_mRNA_fantom_prepare.csv

## step4.6 peak 长度的统计 
cat robust_BED_feature.gff3 | awk '{print $5-$4}' | python -c 'import sys;import pandas as pd ; a = [int(i) for i in sys.stdin];b = pd.Series(a);print(b.describe())' > peak_length.csv
## step4.7 peak 相距的距离统计 
>region_space
for chrom in 1 2 3 4 5 C M
do 
    cat robust_BED_feature.gff3 | grep Chr$chrom | awk '{printf $5","$4","}' | awk -F',' '{for(i=2;i<=(NF-1);i++){print $i}}' | python -c 'import sys;a = [int(i) for i in sys.stdin];c = [print(a[i]-a[i-1]) for i in range(1, len(a),2)]' >> region_space
done
cat region_space | python -c 'import sys;import pandas as pd ; a = [int(i) for i in sys.stdin];b = pd.Series(a);print(b.describe())' > gap_length.csv
## step4.7 TSS_cage_peak_region_stat.R
cat robust_BED_feature.gff3 | awk '{print $5-$4}' | sort |uniq -c | awk '{print $2","$1}' > peaks_region_length_stat.csv
## step4.8 TSS_robust_peak_per_gene.R
bedtools intersect -a robust_BED_feature.gff3 -b $OUTDIR/01.refer_stat/Araport11_GFF3_gene -wb | awk '{print $18}' | awk -F';' '{print $1}' | sed 's/ID=//' | sort | uniq -c | awk 'BEGIN{print "gene,RobustPeakCount"}{print $2","$1}' > robust_peaks_per_gene


mkdir $OUTDIR/05.model_data
cd $OUTDIR/05.model_data
# step5 Robust TSS for model training ----------
less $OUTDIR/01.refer_stat/Araport11_GFF3_mRNA | awk -F "\t" '{OFS="\t";{print $1,$2, $3, $4, $4+1, $6, $7, $8,$9}}' > genome_anotated_tss
# step5.1 tss one point # 2269
bedtools intersect -a $OUTDIR/04.robust/robust_BED_feature.gff3 -b genome_anotated_tss -wb > tss_cover.csv
# step5.2 tss range [-10,10] # 3271
cat $OUTDIR/01.refer_stat/Araport11_GFF3_mRNA | awk -F "\t" '{OFS="\t";if($4-10>0){start=$4-10}else{start=1};{print $1, $2, $3, start, $4+10, $6, $7, $8, $9}}' > genome_anotated_tss_range10
bedtools intersect -a $OUTDIR/04.robust/robust_BED_feature.gff3 -b genome_anotated_tss_range10 -wb > range10_tss_cover.csv
# step5.3 tss range [-50,50] # 7810 
less $OUTDIR/01.refer_stat/Araport11_GFF3_mRNA | awk -F "\t" '{OFS="\t";if($4-50>0){start=$4-50}else{start=1};{print $1, $2, $3, start, $4+50, $6, $7, $8, $9}}' > genome_anotated_tss_range50
bedtools intersect -a $OUTDIR/04.robust/robust_BED_feature.gff3 -b genome_anotated_tss_range50 -wb > range50_tss_cover.csv
# step5.4 tss range [-100,100] # 12911
less $OUTDIR/01.refer_stat/Araport11_GFF3_mRNA | awk -F "\t" '{OFS="\t";if($4-100>0){start=$4-100}else{start=1};{print $1,$2, $3, start, $4+100, $6, $7, $8, $9}}' > genome_anotated_tss_range100
bedtools intersect -a $OUTDIR/04.robust/robust_BED_feature.gff3 -b genome_anotated_tss_range100 -wb > range100_tss_cover.csv
# step5.5 tss range [-500,500] # 42542
less $OUTDIR/01.refer_stat/Araport11_GFF3_mRNA | awk -F "\t" '{OFS="\t";if($4-500>0){start=$4-500}else{start=1};{print $1, $2, $3, start, $4+500, $6, $7, $8, $9}}' > genome_anotated_tss_range500
bedtools intersect -a $OUTDIR/04.robust/robust_BED_feature.gff3 -b genome_anotated_tss_range500 -wb > range500_tss_cover.csv
cat range500_tss_cover.csv | awk '{print $18}' | awk -F';' '{print $1}' | sed 's/ID=//' | sort | uniq -c | awk '{print $2}' > tss_true.csv

# step6 sharp and board promoters
mkdir $OUTDIR/06.sharp_board
cd $OUTDIR/06.sharp_board
cat ${REFER_GFF} | awk '$3=="five_prime_UTR"' > five_prime_UTR.gff
bedtools intersect -a $OUTDIR/04.robust/robust_BED_feature.gff3 -b five_prime_UTR.gff -wb > sharp_board.csv

# step7
mkdir $OUTDIR/07.permissive
cd $OUTDIR/07.permissive
cp $PERMISSIVE .
gunzip tc.spi_merged.ctssMaxCounts3.bed.gz
gt bed_to_gff3 -force -o tc.spi_merged.ctssMaxCounts3.gff3 tc.spi_merged.ctssMaxCounts3.bed

bedtools intersect -a tc.spi_merged.ctssMaxCounts3.gff3 -b ${REFER_GFF} -wb > permissive.bbb
>cage_gff_intersect_permissive.txt
for feature in gene mRNA five_prime_UTR;
do
   echo $feature `cat permissive.bbb | awk '$12=="'$feature'"' | awk '{print $18}' | awk -F';' '{print $1}' | sed 's/ID=//' | sort | uniq | wc -l` >> cage_gff_intersect_permissive.txt
done

# step8 precise five_prime_UTR
mkdir $OUTDIR/08.precise
cd $OUTDIR/08.precise
## peak in utr5 gene numbers
cat $OUTDIR/07.permissive/permissive.bbb | awk '$12=="five_prime_UTR"' | awk '{print $18}' | awk -F';' '{print $1}' | sed 's/ID=//' | sort | uniq | awk -F':' '{print $1}' | sort | uniq > permissive_utr5_genelist


# step9 jbrowser visualization
mkdir $OUTDIR/09.jbrowser
cd $OUTDIR/09.jbrowser
## step9.1 p1@gene_name
bedtools intersect -a $OUTDIR/01.refer_stat/Araport11_GFF3_five_prime_UTR -b $OUTDIR/07.permissive/tc.spi_merged.ctssMaxCounts3.bed -wb > p_format_in_fantom_prepare_permissive.csv

## step9.2 非5' utr区间 p@locus
bedtools intersect -a $OUTDIR/01.refer_stat/Araport11_GFF3_mRNA -b $OUTDIR/07.permissive/tc.spi_merged.ctssMaxCounts3.bed -wb > p_format_in_mRNA_fantom_prepare_permissive.csv
ln -s $OUTDIR/04.robust/p_format_in_mRNA_fantom_prepare.csv  
ln -s $OUTDIR/04.robust/p_format_in_fantom_prepare.csv

python $OUTDIR/../scripts/get_p_format_bed.py
bgzip robust.sort.bed
tabix robust.sort.bed.gz
bgzip permissive.sort.bed
tabix permissive.sort.bed.gz

# step11 在基因上，peak 密度和转录本密度的比较
mkdir $OUTDIR/10.density
cd $OUTDIR/10.density
cat $OUTDIR/01.refer_stat/Araport11_GFF3_mRNA | cut -f9 | awk -F';' '{print $2}' | sed 's/Parent=//g' | sort | uniq -c | awk 'BEGIN{print "Gene,TranscripNum"}{print $2","$1}' > transcripts_per_gene.csv
cat $OUTDIR/09.jbrowser/robust_utr5.bed |cut -f4|awk -F'@' '{print $2}' | sort | uniq -c | awk 'BEGIN{print "Gene,PeakNum"}{print $2","$1}' > peaks_per_gene.csv
python $OUTDIR/../scripts/peak_transcript_density.py
