import subprocess
import glob

# 生成文件列表
sra_files = glob.glob('*.sra')

for sra_file in sra_files:
    base_name = sra_file[:11]
    fq1 = base_name + '_1.fastq'
    fq2 = base_name + '_2.fastq'
    bam_in = base_name + '.sam'
    bam_out = base_name + '.bam'
    sort_out = base_name + '.sorted.bam'
    txt_out = base_name + 'count.txt'

    # 解压SRA文件
    subprocess.run(['fasterq-dump', '--split-3', '-e', '8', sra_file])

    # HISAT2比对
    subprocess.run(['hisat2', '-x', '/home/dell/reference/GRCH38', '-p', '40', '-1', fq1, '-2', fq2, '-S', bam_in])

    # 删除FASTQ文件
    subprocess.run(['rm', '-rf', fq1, fq2], check=True)

    # SAMTOOLS命令
    subprocess.run(['samtools', 'view', '-bS', bam_in, '-o', bam_out])
    subprocess.run(['samtools', 'sort', '-l', '9', '-m', '5G', '-O', 'bam', '-T', 'tmp', '-@', '40', bam_out, '-o', sort_out])
    subprocess.run(['samtools', 'index', sort_out])

    # 删除中间文件
    subprocess.run(['rm', '-rf', bam_in, bam_out], check=True)

    # featureCounts
    subprocess.run(['featureCounts', '-T', '16', '-p', '-t', 'similarity', '-a', '/home/dell/reference/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.out.gtf', '-o', txt_out, sort_out])

sorted_bams = ' '.join(glob.glob('*.sorted.bam'))

# 然后，正确执行featureCounts命令
subprocess.run(f'featureCounts -T 16 -p -t similarity -a /home/dell/reference/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.out.gtf -o all_LTR.txt {sorted_bams}', shell=True)
subprocess.run(f'featureCounts -T 16 -p -a /home/dell/reference/Homo_sapiens.GRCh38.109.chr_patch_hapl_scaff.gtf -o all_exon.txt {sorted_bams}', shell=True)