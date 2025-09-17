#!/bin/sh
#SBATCH -o job.%j.ns.%N.out
#SBATCH --partition=cpu
#SBATCH -J vhost_lt
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=7
#SBATCH --mail-type=end

#####所需文件1.vOTU /lustre/home/liutang/31YellowRiver_meta/06_cdhit_checkv/cdhit/combined.fa 
#####所需文件2.宏基因组组装bins:MAGs /lustre/home/liutang/31YellowRiver_meta/23Binning_YellowRiver_metagenomes
#####共四种方法，每种方法生成一个votu-host 预测结果文件METHOD_votu_mags.txt，展示votu host 关系对
#####所涉及软件，除php运行较慢，其他速度均较快
###为bins做标记
module load blast+/2.9.0-gcc-4.8.5

votu_path=/lustre/home/liutang_faculty/19CBHv/02/22cdhit
binning_path=/lustre/home/liutang_faculty/19CBHv/02/06drep/res_merged/dereplicated_genomes
#binning_path=/lustre/home/liutang_faculty/11basic_lake/03NM/10binning
gtdb=/lustre/home/liutang_faculty/19CBHv/02/07GTDB/res
wp=${PWD}

#####方法3 crispr match 
mkdir 03crispr
#conda activate  /lustre/home/liutang_faculty/10software/mambaforge/ens/R4.2
#minced -spacers 00all_bin_marked.fa 
module load anaconda3
source activate /lustre/home/liutang_faculty/10software/mambaforge/envs/PHP
/lustre/home/liutang_faculty/10software/minced-master/minced -spacers 00all_bin_marked.fa 03crispr/all_bin_marked.crisprs 03crispr/all_bin_marked.gff

makeblastdb -in ${votu_path}/virus_cdhit.fa  -dbtype nucl -out 03crispr/all_votu -hash_index -max_file_sz '2GB'

blastn -db 03crispr/all_votu  \
  -query 03crispr/all_bin_marked_spacers.fa \
  -out 03crispr/all_spacers_votu.blastn.out \
  -outfmt '6 qseqid sseqid pident evalue mismatch bitscore length qlen qcovs qcovhsp slen gapopen qstart qend sstart send' \
  -num_threads ${SLURM_NTASKS} -evalue 1e-5


##按阈值筛选
cat 03crispr/all_spacers_votu.blastn.out |awk '{if($5<=1 && $9==100) print $2,$1}'| awk -F '_[kN]' '{print $1}'| sort|uniq > final03_crispr_votu_mags.txt


#####方法4 Oligonucleotide frequency(ONF) 
#mkdir 04ONF
#cd 04ONF

###软件 php,对于大数据量运行较慢，可以先跑中高质量votu，需要每条votu一个.fa文件 
###将votu分割成1条序列1个文件
#mkdir chm ##新建存放单条序列.fa文件的目录
#cd chm
#awk '/^>/{s=++num}{print > s".fa"}' ${votu_path}/virus_cdhit.fa
#cd ..


#host=${binning_path}/all_bins ##所有bins所在目录
#output=${root_path}/02ONF ##输出路径
#virus=${root_path}/chm ###所有votu单条序列文件所在目录

#source activate PHP
#python /lustre/home/liutang/01software/PHP/PHP-master/countKmer.py -f ${host} -d ./ -n all_bins_hostkmer -c 1

#python /lustre/home/liutang/01software/PHP/PHP-master/PHP.py -v chm -o ./ -d ./ -n all_bins_hostkmer

#将单序列.fa文件与votu id 对应起来
#for i in  chm/*
#do
#name=${i##*/}
#sq=$(awk -F '>' 'NR==1{print $2}' $i)
#echo $name,$sq >> num_sq.txt
#done
#得到votu-host文件，选用阈值score>1436
#cat all_bins_hostkmer_Prediction_Maxhost.tsv |awk '{if($2>1436 && NR >1 )print $0}' > max_1436.tsv
#python /lustre/home/liutang/31Yellowriver_virus/18virhost/10php/yellowriver/findhost.py ./ ##替换output路径
#cat 1436_prediction.txt | awk '{print $1,$3}'|sort|uniq > ${root_path}/onf_votu_mags.txt



#####合并所有方法结果
#cat ss_votu_mags.txt onf_votu_mags.txt trna_votu_mags.txt crispr_votu_mags.txt |sort|uniq > 4meth_votu_mags.txt
cat final01_votu_mags.txt  final02_rna_votu_mags.txt  final03_crispr_votu_mags.txt |sort|uniq > all_votu_mags.txt
#####添加bins的物种注释信息
python /lustre/home/liutang_faculty/03code/lt/addclas.py ./ ##替换4meth_votu_mags.txt所在路径
cat final_pair_clas.txt | awk 'BEGIN{FS=OFS="\t"}{print $1FS$3}' | sort|uniq > final_vir_clas_uniq.txt
#####统计跨门感染情况
cat final_vir_clas_uniq.txt |awk -F ';' '{print $1,$2}'|sort|uniq  > final_votu_phylum.txt
cat final_votu_phylum.txt | awk '{print $1}' |sort|uniq -c  |awk '{if($1>1) print $2}' > final_votu_multi_phylum.txt





