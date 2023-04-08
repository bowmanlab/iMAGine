#!/bin/bash

# author: Avishek Dutta, avdutta@ucsd.edu

helpFunction()
{
   echo -e "\t Usage: ./imagine.sh -r LMG2001_200121_7-R1.fastq.gz -R LMG2001_200121_7-R2.fastq.gz -s sample_name |& tee -a imagine.txt"
   echo -e "\t-r read 1"
   echo -e "\t-R read 2"
   echo -e "\t-s name of the sample (modifier) given to intermediated files"
   echo -e "\t-h help"
   echo -e "Modify MEMORY and CPU variables in imagine.sh as appropriate"
   exit 1 # Exit script after printing help
}

while getopts r:R:s:h opt
do
   case "${opt}" in
      r ) read1=${OPTARG};;
      R ) read2=${OPTARG};;
      s ) name=${OPTARG};;
      h ) helpFunction ;;
   esac
done

START=$(date +%s)

### testing ###

#read1=/data_store/seq_data/ucsd_microbiome_20230214/Aluwihare_DNA6_S10_L001_R1_001.fastq.gz
#read2=/data_store/seq_data/ucsd_microbiome_20230214/Aluwihare_DNA6_S10_L001_R2_001.fastq.gz
#name=Aluwihare_DNA6

### running ###

##

MEMORY=1000
CPUS=64

fastp -i ${read1} -I ${read2} -o ${name}_filt_R1.fq.gz -O ${name}_filt_R2.fq.gz -e 30 -w 1 -j fastp_${name}.json -h fastp_${name}.html &&

## --only-assembler flag used because error correction seems to hang for many samples.  The --continue line is provided
## here in case you need to recover from a stop.

metaspades.py -1 ${name}_filt_R1.fq.gz -2 ${name}_filt_R2.fq.gz -k 21,33,55 -o metaspades_output_${name} -t ${CPUS} -m ${MEMORY} --only-assembler &&

#metaspades.py -1 ${name}_filt_R1.fq.gz -2 ${name}_filt_R2.fq.gz -k 21,33,55 -o metaspades_output_${name} -t 64 -m ${MEMORY} &&

#metaspades.py -o metaspades_output_${name} --continue &&

quast.py metaspades_output_${name}/contigs.fasta -o quast_output_${name} &&

mkdir binning_${name} &&

cp metaspades_output_${name}/contigs.fasta binning_${name} &&

cd binning_${name} &&

bwa index contigs.fasta &&

bwa mem -t ${CPUS} contigs.fasta ../${name}_filt_R1.fq.gz ../${name}_filt_R2.fq.gz > ${name}_map.sam &&

samtools view -S -b ${name}_map.sam > ${name}_map.bam &&

rm ${name}_map.sam

samtools sort -o ${name}_map_sorted.bam -O bam ${name}_map.bam &&

jgi_summarize_bam_contig_depths --outputDepth ${name}_depth.txt ${name}_map_sorted.bam &&

metabat2 -i contigs.fasta -a ${name}_depth.txt -o bins_dir/${name}_bin --seed 1234 -m 1500 -v &&

checkm data setRoot /home/jsbowman/checkm_files &&

checkm lineage_wf bins_dir/ checkm/ -x .fa -t ${CPUS} &&

mag_extract.py &&

export GTDBTK_DATA_PATH=/home/jsbowman/gtdbtk_files &&

gtdbtk identify --genome_dir high_qual_draft/ --out_dir high_qual_identify_output --cpus ${CPUS} --extension fa
gtdbtk identify --genome_dir medium_qual_draft/ --out_dir medium_qual_identify_output --cpus ${CPUS} --extension fa
gtdbtk identify --genome_dir low_qual_draft/ --out_dir low_qual_identify_output --cpus ${CPUS} --extension fa

gtdbtk align --identify_dir high_qual_identify_output --out_dir high_qual_align_output --cpus ${CPUS}
gtdbtk align --identify_dir medium_qual_identify_output --out_dir medium_qual_align_output --cpus ${CPUS}
gtdbtk align --identify_dir low_qual_identify_output --out_dir low_qual_align_output --cpus ${CPUS}

gtdbtk classify --genome_dir high_qual_draft/ --align_dir high_qual_align_output/ --out_dir high_qual_classify_output --cpus ${CPUS} --extension fa
gtdbtk classify --genome_dir medium_qual_draft/ --align_dir medium_qual_align_output/ --out_dir medium_qual_classify_output --cpus ${CPUS} --extension fa
gtdbtk classify --genome_dir low_qual_draft/ --align_dir low_qual_align_output/ --out_dir low_qual_classify_output --cpus ${CPUS} --extension fa

gtdbtk ani_rep --genome_dir high_qual_draft/ --out_dir high_qual_ani_rep/ --cpus ${CPUS} --extension fa
gtdbtk ani_rep --genome_dir medium_qual_draft/ --out_dir medium_qual_ani_rep/ --cpus ${CPUS} --extension fa
gtdbtk ani_rep --genome_dir low_qual_draft/ --out_dir low_qual_ani_rep/ --cpus ${CPUS} --extension fa

cd ../

END=$(date +%s)
DIFF=$(( $END - $START ))

echo "It took $DIFF seconds"
