import os
from tqdm import *
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--path')
args = parser.parse_args()
ref_path = "/home/undergraduate/reference/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa"
original_path = args.path
folders = os.listdir(original_path)

for i in tqdm(folders):
  try:
    fq = os.listdir(os.path.join(original_path,i))
    fq = [n for n in fq if n.split('.')[-1] != 'txt']
    fq1 = original_path+'/'+i+'/'+fq[0]
    fq2 = original_path+'/'+i+'/'+fq[1]
    outpath = original_path + '/' + i
    os.system(f'trim_galore -j 1 -q 20 --phred33 --stringency 3 --length 20 -e 0.1 --paired {fq1} {fq2} -o {outpath} --gzip --cores 4')
    fq = os.listdir(os.path.join(original_path,i))
    fq = [n for n in fq if '_val_' in n]
    fq1 = original_path+'/'+i+'/'+fq[0]
    fq2 = original_path+'/'+i+'/'+fq[1]
    sampath = original_path + '/' + i + '/out.sam'
    bampath = original_path + '/' + i + '/out.bam'
    sortpath = original_path + '/' + i + '/out.sort.bam'
    duppath = original_path + '/' + i + '/out.dup.bam'
    rawpath = original_path + '/' + i + '/raw.vcf'
    filtered = original_path + '/' + i + '/filtered.vcf'
    anno = original_path + '/' + i + '/filtered.anno.vcf'
    os.system(f'bwa-mem2 mem -t 60  {ref_path} {fq1} {fq2} > {sampath}')
    os.system(f'samtools view -@ 15 -bf 3 {sampath} > {bampath} ')
    os.system(f'samtools sort -@ 15 {bampath} > {sortpath}')
    os.system(f'java -jar ~/picard.jar MarkDuplicates --REMOVE_DUPLICATES true --OPTICAL_DUPLICATE_PIXEL_DISTANCE 100 -I {sortpath} -M dup.txt -O {duppath}')
    os.system(f'freebayes -f {ref_path} {duppath} > {rawpath}')
    os.system(f'vcffilter -f "QUAL > 20 & DP > 5 & AF > 0.8" {rawpath} > {filtered}')
    os.system(f'java -jar ~/software/snpEff/snpEff.jar eff WBcel235.99 {filtered} > {anno}')
    os.system('rm out.sam')
    os.system('rm out.sort.bam')
  except:
    print(i)