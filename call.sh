#!/bin/sh

# $1, $2 - fastq files
genomeDir=/data9/bio/runs-ms/human_ref/hg19/STAR

dir=/data9/bio/TCGA_data/COAD_RNA/FFPE_calling/frozen/$(basename $1 _primary_tumor_R1.fastq.gz)
mkdir -p $dir
runDir=$dir/STAR_1pass

# mapping 1st pass
mkdir -p $runDir
cd $runDir
echo -e "Started STAR 1st pass at $(date)\n"
STAR --genomeDir $genomeDir --readFilesIn $1 $2 --runThreadN 20 --readFilesCommand zcat
echo -e "Finished STAR 1st pass at $(date)\n"

# creating index for 2nd pass

genomeDir2=$runDir/index
mkdir -p $genomeDir2
echo -e "Started index generation for 2nd pass at $(date)\n"
STAR --runMode genomeGenerate --genomeDir $genomeDir2 --genomeFastaFiles $genomeDir/../hg19.fa --sjdbFileChrStartEnd $runDir/SJ.out.tab --sjdbOverhang 75 --runThreadN 20
echo -e "Finished index generation for 2nd pass at $(date)\n"

# mapping 2nd pass

runDir2=$dir/STAR_2pass
mkdir -p $runDir2
cd $runDir2
echo -e "Started STAR 2nd pass at $(date)\n"
STAR --genomeDir $genomeDir2 --readFilesIn $1 $2 --runThreadN 20 --readFilesCommand zcat
echo -e "Finished STAR 2nd pass at $(date)\n"

# Add read groups, sort, mark duplicates, and create index - DONE
cd $runDir2
echo -e "Started AddOrReplaceReadGroups at $(date)\n"
/srv/common/opt/jdk1.8.0_91/bin/java -jar /data7a/bio/human_genomics/shared/tools/bcbio/anaconda/share/picard-2.18.2-0/picard.jar AddOrReplaceReadGroups I=Aligned.out.sam O=rg_added_sorted.bam SO=coordinate RGID=id RGLB=library RGPL=platform RGPU=machine RGSM=sample
echo -e "Started MarkDuplicates at $(date)\n"
/srv/common/opt/jdk1.8.0_91/bin/java -jar /data7a/bio/human_genomics/shared/tools/bcbio/anaconda/share/picard-2.18.2-0/picard.jar MarkDuplicates I=rg_added_sorted.bam O=dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics

#  Split'N'Trim and reassign mapping qualities - DONE

cd $runDir2
echo -e "Started Split'N'Trim and reassign mapping qualities at $(date)\n"
/srv/common/opt/jdk1.8.0_91/bin/java -jar /data7a/bio/human_genomics/shared/tools/bcbio/anaconda/share/picard-2.18.2-0/picard.jar ReorderSam INPUT=dedupped.bam OUTPUT=karyo.bam REFERENCE=/data9/bio/runs-ms/human_ref/hg19/karyo/ucsc.hg19.fasta
samtools index karyo.bam
/srv/common/opt/jdk1.8.0_91/bin/java -jar /data7a/bio/human_genomics/shared/tools/GATK/3.5-0-g36282e4/GenomeAnalysisTK.jar -T SplitNCigarReads -R $genomeDir/../hg19.fa -I karyo.bam -o split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

echo -e "BAM generation finished at $(date)\n"


echo -e "HaplotypeCaller started at $(date)\n"
mkdir -p $runDir2/hcaller
while read chr
do
        echo "/srv/common/opt/jdk1.8.0_91/bin/java -jar /data7a/bio/human_genomics/shared/tools/GATK/3.5-0-g36282e4/GenomeAnalysisTK.jar -T HaplotypeCaller -R /data9/bio/runs-ms/human_ref/hg19/karyo/ucsc.hg19.fasta -I $runDir2/split.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -L $chr -o $runDir2/hcaller/${chr}.vcf" | qsub -cwd -N ${chr}_${name}_hcaller -o qsub_logs/${name}_hcaller_${chr}.out -e qsub_logs/${name}_hcaller_${chr}.err

#/srv/common/opt/jdk1.8.0_91/bin/java -jar /data7a/bio/human_genomics/shared/tools/GATK/3.5-0-g36282e4/GenomeAnalysisTK.jar -T VariantFiltration -R hg_19.fasta -V input.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o output.vcf

done < /data9/bio/runs-ms/human_ref/hg19/karyo/chr_names

