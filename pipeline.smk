configfile: "./9.config.yaml"
workdir: "/media/nvme8T/WGS_pipeline_storage"

rule all:
    input:
        expand('/media/nvme8T/WGS_pipeline_storage/{sample}.VQSR.vcf.gz',sample = config["samples"]),
#        expand('/media/nvme8T/WGS_pipeline_storage/{sample}_qualimap_result/{sample}.pdf',sample = config["samples"]),
        expand('/media/nvme8T/WGS_pipeline_storage/{sample}.BQSR.soft_clip.hg38noalt.bam',sample = config["samples"])
#        expand('{sample}.g.chrX.vcf.gz',sample = config["samples"])

def getname(wildcards):
    return config["samples"][wildcards.sample]

rule fastp:
    input:
        r1= '/media/nvme8T/WGS_pipeline/{sample}/{sample}.R1.fastq.gz',
        r2= '/media/nvme8T/WGS_pipeline/{sample}/{sample}.R2.fastq.gz'
    output:
        r1= temp('/media/nvme8T/WGS_pipeline/{sample}/{sample}.fastp.R1.fastq.gz'),
        r2= temp('/media/nvme8T/WGS_pipeline/{sample}/{sample}.fastp.R2.fastq.gz')
    params:
        htmlname = '/media/nvme8T/WGS_pipeline_storage/{sample}.html',
        jsonname = '/media/nvme8T/WGS_pipeline_storage/{sample}.json'
    threads:30
    log:
        '/media/nvme8T/WGS_pipeline/{sample}/{sample}.fastp.log'
    shell:
        '(/home/xs/genome_tool/fastp/fastp.0.23.2 -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} -h {params.htmlname} -j {params.jsonname} -w {threads}) 2>{log}'        

rule bwa_mem2_to_bam:
    input:
        r1= '/media/nvme8T/WGS_pipeline/{sample}/{sample}.fastp.R1.fastq.gz',
        r2= '/media/nvme8T/WGS_pipeline/{sample}/{sample}.fastp.R2.fastq.gz'
    output:
        temp('/media/nvme8T/WGS_pipeline/{sample}/{sample}.bam')
    params:
        name = getname,
        ref = '/media/nvme1/reference/hg38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna',
#        rg = "'@RG\\tID:{params.name}\\tPL:Illumina\\tPU:{params.name}\\tLB:{params.name}\\tSM:{params.name}'"
    threads: 46
    log:
        '/media/nvme8T/WGS_pipeline/{sample}/{sample}.bwa.log'
    resources: tmpdir='/media/nvme8T/temp'
    shell:
        '(/home/xs/genome_tool/bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem -Y -t {threads} -R '
        '\'@RG\\tID:{params.name}\\tPL:Illumina\\tPU:{params.name}\\tLB:{params.name}\\tSM:{params.name}\' {params.ref} '
        '{input.r1} {input.r2} | samtools view -S -b - > {output}) 2>{log}'

rule sort_dup_mark:
    input:
        '/media/nvme8T/WGS_pipeline/{sample}/{sample}.bam'
    output:
        r1 = temp('/media/nvme8T/WGS_pipeline/{sample}/{sample}.dup.bam'),
        r2 = temp('/media/nvme8T/WGS_pipeline/{sample}/{sample}.dup.bam.bai'),
        r3 = temp('/media/nvme8T/WGS_pipeline/{sample}/{sample}.dup.bam.sbi')
    threads: 50
    params:
        name = getname
    log:
        '/media/nvme8T/WGS_pipeline/{sample}/{sample}.sort_dup_mark.log'
    resources: tmpdir='/media/nvme8T/temp'
    shell:
        '(/home/xs/genome_tool/gatk-4.6.2.0/gatk MarkDuplicatesSpark -I {input} -O {output.r1} '
        '--spark-master local[{threads}] --conf "spark.local.dir=/media/xs/4t_cache1/spark_temp") 2> {log}'

rule settags:
    input:
        r1 = '/media/nvme8T/WGS_pipeline/{sample}/{sample}.dup.bam',
        r2 = '/media/nvme8T/WGS_pipeline/{sample}/{sample}.dup.bam.bai',
        r3 = '/media/nvme8T/WGS_pipeline/{sample}/{sample}.dup.bam.sbi'
    output:
        temp('/media/nvme8T/WGS_pipeline/{sample}/{sample}.dup.fix.bam')
    params:
        ref = '/media/nvme1/reference/hg38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna'
    log:
        '/media/nvme8T/WGS_pipeline/{sample}/{sample}.settags.log'
    resources: tmpdir='/media/nvme8T/temp'
    shell:
        '(/home/xs/genome_tool/gatk-4.6.2.0/gatk SetNmMdAndUqTags '
        '-R {params.ref} -I {input.r1} -O {output} --TMP_DIR /media/nvme8T/temp) 2>{log}'

rule BQSR:
    input:
        '/media/nvme8T/WGS_pipeline/{sample}/{sample}.dup.fix.bam'
    output:
        '/media/nvme8T/WGS_pipeline_storage/{sample}.BQSR.soft_clip.hg38noalt.bam'
    params:
        ref = '/media/nvme1/reference/hg38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna',
        SNPS_1000G = '/home/xs/genome_tool/reference-genome/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz',
        indels1 = '/home/xs/genome_tool/reference-genome/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz',
        indels2 = '/home/xs/genome_tool/reference-genome/hg38/Homo_sapiens_assembly38.known_indels.vcf.gz'
    log:
        '/media/nvme8T/WGS_pipeline/{sample}/{sample}.BQSR.log'
    threads: 50
    resources: tmpdir='/media/nvme8T/temp'
    shell:
        '(/home/xs/genome_tool/gatk-4.6.2.0/gatk BQSRPipelineSpark --java-options "-Dsamjdk.compression_level=5" '
        '-R {params.ref} '
        '--known-sites {params.SNPS_1000G} --known-sites {params.indels1} --known-sites {params.indels2} '
        '-I {input} -O {output} --spark-master local[{threads}] '
        '--conf "spark.local.dir=/media/xs/4t_cache1/spark_temp") 2>{log}'



rule HaplotypeCaller_chrY:
    input:
        r1 = '/media/nvme8T/WGS_pipeline_storage/{sample}.BQSR.soft_clip.hg38noalt.bam'
    output:
        temp('/media/nvme8T/WGS_pipeline/{sample}/{sample}.g.chrY.vcf.gz')
    params:
        ref = '/media/nvme1/reference/hg38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna',
        dbsnp = '/home/xs/genome_tool/reference-genome/hg38/dbsnp_146.hg38.vcf.gz'
    log:
        '/media/nvme8T/WGS_pipeline/{sample}/{sample}.chrY.HaplotypeCaller.log'
    resources: tmpdir='/media/nvme8T/temp'
    threads: 2
    shell:
        '(/home/xs/genome_tool/gatk-4.6.2.0/gatk HaplotypeCaller -R {params.ref} '
        '--native-pair-hmm-threads {threads} --intervals chrY '
        '-I {input.r1} -D {params.dbsnp} '
        '-O {output} --tmp-dir /media/nvme8T/temp --ERC GVCF) 2>{log}'

rule HaplotypeCaller_chrX:
    input:
        '/media/nvme8T/WGS_pipeline_storage/{sample}.BQSR.soft_clip.hg38noalt.bam'
    output:
        temp('/media/nvme8T/WGS_pipeline/{sample}/{sample}.g.chrX.vcf.gz')
    params:
        ref = '/media/nvme1/reference/hg38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna',
        dbsnp = '/home/xs/genome_tool/reference-genome/hg38/dbsnp_146.hg38.vcf.gz'
    log:
        '/media/nvme8T/WGS_pipeline/{sample}/{sample}.chrX.HaplotypeCaller.log'
    resources: tmpdir='/media/nvme8T/temp'
    threads: 2
    shell:
        '(/home/xs/genome_tool/gatk-4.6.2.0/gatk HaplotypeCaller -R {params.ref} '
        '--native-pair-hmm-threads {threads} --intervals chrX '
        '-I {input} -D {params.dbsnp} '
        '-O {output} --tmp-dir /media/nvme8T/temp --ERC GVCF) 2>{log}'

rule HaplotypeCaller_chr1:
    input:
        '/media/nvme8T/WGS_pipeline_storage/{sample}.BQSR.soft_clip.hg38noalt.bam'
    output:
        temp('/media/nvme8T/WGS_pipeline/{sample}/{sample}.g.chr1.vcf.gz')
    params:
        ref = '/media/nvme1/reference/hg38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna',
        dbsnp = '/home/xs/genome_tool/reference-genome/hg38/dbsnp_146.hg38.vcf.gz'
    log:
        '/media/nvme8T/WGS_pipeline/{sample}/{sample}.chr1.HaplotypeCaller.log'
    resources: tmpdir='/media/nvme8T/temp'
    threads: 2
    shell:
        '(/home/xs/genome_tool/gatk-4.6.2.0/gatk HaplotypeCaller -R {params.ref} '
        '--native-pair-hmm-threads {threads} --intervals chr1 '
        '-I {input} -D {params.dbsnp} '
        '-O {output} --tmp-dir /media/nvme8T/temp --ERC GVCF) 2>{log}'

rule HaplotypeCaller_chr2:
    input:
        '/media/nvme8T/WGS_pipeline_storage/{sample}.BQSR.soft_clip.hg38noalt.bam'
    output:
        temp('/media/nvme8T/WGS_pipeline/{sample}/{sample}.g.chr2.vcf.gz')
    params:
        ref = '/media/nvme1/reference/hg38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna',
        dbsnp = '/home/xs/genome_tool/reference-genome/hg38/dbsnp_146.hg38.vcf.gz'
    log:
        '/media/nvme8T/WGS_pipeline/{sample}/{sample}.chr2.HaplotypeCaller.log'
    resources: tmpdir='/media/nvme8T/temp'
    threads: 2
    shell:
        '(/home/xs/genome_tool/gatk-4.6.2.0/gatk HaplotypeCaller -R {params.ref} '
        '--native-pair-hmm-threads {threads} --intervals chr2 '
        '-I {input} -D {params.dbsnp} '
        '-O {output} --tmp-dir /media/nvme8T/temp --ERC GVCF) 2>{log}'

rule HaplotypeCaller_chr3:
    input:
        '/media/nvme8T/WGS_pipeline_storage/{sample}.BQSR.soft_clip.hg38noalt.bam'
    output:
        temp('/media/nvme8T/WGS_pipeline/{sample}/{sample}.g.chr3.vcf.gz')
    params:
        ref = '/media/nvme1/reference/hg38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna',
        dbsnp = '/home/xs/genome_tool/reference-genome/hg38/dbsnp_146.hg38.vcf.gz'
    log:
        '/media/nvme8T/WGS_pipeline/{sample}/{sample}.chr3.HaplotypeCaller.log'
    resources: tmpdir='/media/nvme8T/temp'
    threads: 2
    shell:
        '(/home/xs/genome_tool/gatk-4.6.2.0/gatk HaplotypeCaller -R {params.ref} '
        '--native-pair-hmm-threads {threads} --intervals chr3 '
        '-I {input} -D {params.dbsnp} '
        '-O {output} --tmp-dir /media/nvme8T/temp --ERC GVCF) 2>{log}'

rule HaplotypeCaller_chr4:
    input:
        '/media/nvme8T/WGS_pipeline_storage/{sample}.BQSR.soft_clip.hg38noalt.bam'
    output:
        temp('/media/nvme8T/WGS_pipeline/{sample}/{sample}.g.chr4.vcf.gz')
    params:
        ref = '/media/nvme1/reference/hg38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna',
        dbsnp = '/home/xs/genome_tool/reference-genome/hg38/dbsnp_146.hg38.vcf.gz'
    log:
        '/media/nvme8T/WGS_pipeline/{sample}/{sample}.chr4.HaplotypeCaller.log'
    resources: tmpdir='/media/nvme8T/temp'
    threads: 2
    shell:
        '(/home/xs/genome_tool/gatk-4.6.2.0/gatk HaplotypeCaller -R {params.ref} '
        '--native-pair-hmm-threads {threads} --intervals chr4 '
        '-I {input} -D {params.dbsnp} '
        '-O {output} --tmp-dir /media/nvme8T/temp --ERC GVCF) 2>{log}'

rule HaplotypeCaller_chr5:
    input:
        '/media/nvme8T/WGS_pipeline_storage/{sample}.BQSR.soft_clip.hg38noalt.bam'
    output:
        temp('/media/nvme8T/WGS_pipeline/{sample}/{sample}.g.chr5.vcf.gz')
    params:
        ref = '/media/nvme1/reference/hg38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna',
        dbsnp = '/home/xs/genome_tool/reference-genome/hg38/dbsnp_146.hg38.vcf.gz'
    log:
        '/media/nvme8T/WGS_pipeline/{sample}/{sample}.chr5.HaplotypeCaller.log'
    resources: tmpdir='/media/nvme8T/temp'
    threads: 2
    shell:
        '(/home/xs/genome_tool/gatk-4.6.2.0/gatk HaplotypeCaller -R {params.ref} '
        '--native-pair-hmm-threads {threads} --intervals chr5 '
        '-I {input} -D {params.dbsnp} '
        '-O {output} --tmp-dir /media/nvme8T/temp --ERC GVCF) 2>{log}'

rule HaplotypeCaller_chr6:
    input:
        '/media/nvme8T/WGS_pipeline_storage/{sample}.BQSR.soft_clip.hg38noalt.bam'
    output:
        temp('/media/nvme8T/WGS_pipeline/{sample}/{sample}.g.chr6.vcf.gz')
    params:
        ref = '/media/nvme1/reference/hg38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna',
        dbsnp = '/home/xs/genome_tool/reference-genome/hg38/dbsnp_146.hg38.vcf.gz'
    log:
        '/media/nvme8T/WGS_pipeline/{sample}/{sample}.chr6.HaplotypeCaller.log'
    resources: tmpdir='/media/nvme8T/temp'
    threads: 2
    shell:
        '(/home/xs/genome_tool/gatk-4.6.2.0/gatk HaplotypeCaller -R {params.ref} '
        '--native-pair-hmm-threads {threads} --intervals chr6 '
        '-I {input} -D {params.dbsnp} '
        '-O {output} --tmp-dir /media/nvme8T/temp --ERC GVCF) 2>{log}'

rule HaplotypeCaller_chr7:
    input:
        '/media/nvme8T/WGS_pipeline_storage/{sample}.BQSR.soft_clip.hg38noalt.bam'
    output:
        temp('/media/nvme8T/WGS_pipeline/{sample}/{sample}.g.chr7.vcf.gz')
    params:
        ref = '/media/nvme1/reference/hg38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna',
        dbsnp = '/home/xs/genome_tool/reference-genome/hg38/dbsnp_146.hg38.vcf.gz'
    log:
        '/media/nvme8T/WGS_pipeline/{sample}/{sample}.chr7.HaplotypeCaller.log'
    resources: tmpdir='/media/nvme8T/temp'
    threads: 2
    shell:
        '(/home/xs/genome_tool/gatk-4.6.2.0/gatk HaplotypeCaller -R {params.ref} '
        '--native-pair-hmm-threads {threads} --intervals chr7 '
        '-I {input} -D {params.dbsnp} '
        '-O {output} --tmp-dir /media/nvme8T/temp --ERC GVCF) 2>{log}'

rule HaplotypeCaller_chr8:
    input:
        '/media/nvme8T/WGS_pipeline_storage/{sample}.BQSR.soft_clip.hg38noalt.bam'
    output:
        temp('/media/nvme8T/WGS_pipeline/{sample}/{sample}.g.chr8.vcf.gz')
    params:
        ref = '/media/nvme1/reference/hg38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna',
        dbsnp = '/home/xs/genome_tool/reference-genome/hg38/dbsnp_146.hg38.vcf.gz'
    log:
        '/media/nvme8T/WGS_pipeline/{sample}/{sample}.chr8.HaplotypeCaller.log'
    resources: tmpdir='/media/nvme8T/temp'
    threads: 2
    shell:
        '(/home/xs/genome_tool/gatk-4.6.2.0/gatk HaplotypeCaller -R {params.ref} '
        '--native-pair-hmm-threads {threads} --intervals chr8 '
        '-I {input} -D {params.dbsnp} '
        '-O {output} --tmp-dir /media/nvme8T/temp --ERC GVCF) 2>{log}'

rule HaplotypeCaller_chr9:
    input:
        '/media/nvme8T/WGS_pipeline_storage/{sample}.BQSR.soft_clip.hg38noalt.bam'
    output:
        temp('/media/nvme8T/WGS_pipeline/{sample}/{sample}.g.chr9.vcf.gz')
    params:
        ref = '/media/nvme1/reference/hg38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna',
        dbsnp = '/home/xs/genome_tool/reference-genome/hg38/dbsnp_146.hg38.vcf.gz'
    log:
        '/media/nvme8T/WGS_pipeline/{sample}/{sample}.chr9.HaplotypeCaller.log'
    resources: tmpdir='/media/nvme8T/temp'
    threads: 2
    shell:
        '(/home/xs/genome_tool/gatk-4.6.2.0/gatk HaplotypeCaller -R {params.ref} '
        '--native-pair-hmm-threads {threads} --intervals chr9 '
        '-I {input} -D {params.dbsnp} '
        '-O {output} --tmp-dir /media/nvme8T/temp --ERC GVCF) 2>{log}'

rule HaplotypeCaller_chr10:
    input:
        '/media/nvme8T/WGS_pipeline_storage/{sample}.BQSR.soft_clip.hg38noalt.bam'
    output:
        temp('/media/nvme8T/WGS_pipeline/{sample}/{sample}.g.chr10.vcf.gz')
    params:
        ref = '/media/nvme1/reference/hg38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna',
        dbsnp = '/home/xs/genome_tool/reference-genome/hg38/dbsnp_146.hg38.vcf.gz'
    log:
        '/media/nvme8T/WGS_pipeline/{sample}/{sample}.chr10.HaplotypeCaller.log'
    resources: tmpdir='/media/nvme8T/temp'
    threads: 2
    shell:
        '(/home/xs/genome_tool/gatk-4.6.2.0/gatk HaplotypeCaller -R {params.ref} '
        '--native-pair-hmm-threads {threads} --intervals chr10 '
        '-I {input} -D {params.dbsnp} '
        '-O {output} --tmp-dir /media/nvme8T/temp --ERC GVCF) 2>{log}'

rule HaplotypeCaller_chr11:
    input:
        '/media/nvme8T/WGS_pipeline_storage/{sample}.BQSR.soft_clip.hg38noalt.bam'
    output:
        temp('/media/nvme8T/WGS_pipeline/{sample}/{sample}.g.chr11.vcf.gz')
    params:
        ref = '/media/nvme1/reference/hg38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna',
        dbsnp = '/home/xs/genome_tool/reference-genome/hg38/dbsnp_146.hg38.vcf.gz'
    log:
        '/media/nvme8T/WGS_pipeline/{sample}/{sample}.chr11.HaplotypeCaller.log'
    resources: tmpdir='/media/nvme8T/temp'
    threads: 2
    shell:
        '(/home/xs/genome_tool/gatk-4.6.2.0/gatk HaplotypeCaller -R {params.ref} '
        '--native-pair-hmm-threads {threads} --intervals chr11 '
        '-I {input} -D {params.dbsnp} '
        '-O {output} --tmp-dir /media/nvme8T/temp --ERC GVCF) 2>{log}'

rule HaplotypeCaller_chr12:
    input:
        '/media/nvme8T/WGS_pipeline_storage/{sample}.BQSR.soft_clip.hg38noalt.bam'
    output:
        temp('/media/nvme8T/WGS_pipeline/{sample}/{sample}.g.chr12.vcf.gz')
    params:
        ref = '/media/nvme1/reference/hg38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna',
        dbsnp = '/home/xs/genome_tool/reference-genome/hg38/dbsnp_146.hg38.vcf.gz'
    log:
        '/media/nvme8T/WGS_pipeline/{sample}/{sample}.chr12.HaplotypeCaller.log'
    resources: tmpdir='/media/nvme8T/temp'
    threads: 2
    shell:
        '(/home/xs/genome_tool/gatk-4.6.2.0/gatk HaplotypeCaller -R {params.ref} '
        '--native-pair-hmm-threads {threads} --intervals chr12 '
        '-I {input} -D {params.dbsnp} '
        '-O {output} --tmp-dir /media/nvme8T/temp --ERC GVCF) 2>{log}'

rule HaplotypeCaller_chr13:
    input:
        '/media/nvme8T/WGS_pipeline_storage/{sample}.BQSR.soft_clip.hg38noalt.bam'
    output:
        temp('/media/nvme8T/WGS_pipeline/{sample}/{sample}.g.chr13.vcf.gz')
    params:
        ref = '/media/nvme1/reference/hg38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna',
        dbsnp = '/home/xs/genome_tool/reference-genome/hg38/dbsnp_146.hg38.vcf.gz'
    log:
        '/media/nvme8T/WGS_pipeline/{sample}/{sample}.chr13.HaplotypeCaller.log'
    resources: tmpdir='/media/nvme8T/temp'
    threads: 2
    shell:
        '(/home/xs/genome_tool/gatk-4.6.2.0/gatk HaplotypeCaller -R {params.ref} '
        '--native-pair-hmm-threads {threads} --intervals chr13 '
        '-I {input} -D {params.dbsnp} '
        '-O {output} --tmp-dir /media/nvme8T/temp --ERC GVCF) 2>{log}'

rule HaplotypeCaller_chr14:
    input:
        '/media/nvme8T/WGS_pipeline_storage/{sample}.BQSR.soft_clip.hg38noalt.bam'
    output:
        temp('/media/nvme8T/WGS_pipeline/{sample}/{sample}.g.chr14.vcf.gz')
    params:
        ref = '/media/nvme1/reference/hg38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna',
        dbsnp = '/home/xs/genome_tool/reference-genome/hg38/dbsnp_146.hg38.vcf.gz'
    log:
        '/media/nvme8T/WGS_pipeline/{sample}/{sample}.chr14.HaplotypeCaller.log'
    resources: tmpdir='/media/nvme8T/temp'
    threads: 2
    shell:
        '(/home/xs/genome_tool/gatk-4.6.2.0/gatk HaplotypeCaller -R {params.ref} '
        '--native-pair-hmm-threads {threads} --intervals chr14 '
        '-I {input} -D {params.dbsnp} '
        '-O {output} --tmp-dir /media/nvme8T/temp --ERC GVCF) 2>{log}'

rule HaplotypeCaller_chr15:
    input:
        '/media/nvme8T/WGS_pipeline_storage/{sample}.BQSR.soft_clip.hg38noalt.bam'
    output:
        temp('/media/nvme8T/WGS_pipeline/{sample}/{sample}.g.chr15.vcf.gz')
    params:
        ref = '/media/nvme1/reference/hg38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna',
        dbsnp = '/home/xs/genome_tool/reference-genome/hg38/dbsnp_146.hg38.vcf.gz'
    log:
        '/media/nvme8T/WGS_pipeline/{sample}/{sample}.chr15.HaplotypeCaller.log'
    resources: tmpdir='/media/nvme8T/temp'
    threads: 2
    shell:
        '(/home/xs/genome_tool/gatk-4.6.2.0/gatk HaplotypeCaller -R {params.ref} '
        '--native-pair-hmm-threads {threads} --intervals chr15 '
        '-I {input} -D {params.dbsnp} '
        '-O {output} --tmp-dir /media/nvme8T/temp --ERC GVCF) 2>{log}'

rule HaplotypeCaller_chr16:
    input:
        '/media/nvme8T/WGS_pipeline_storage/{sample}.BQSR.soft_clip.hg38noalt.bam'
    output:
        temp('/media/nvme8T/WGS_pipeline/{sample}/{sample}.g.chr16.vcf.gz')
    params:
        ref = '/media/nvme1/reference/hg38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna',
        dbsnp = '/home/xs/genome_tool/reference-genome/hg38/dbsnp_146.hg38.vcf.gz'
    log:
        '/media/nvme8T/WGS_pipeline/{sample}/{sample}.chr16.HaplotypeCaller.log'
    resources: tmpdir='/media/nvme8T/temp'
    threads: 2
    shell:
        '(/home/xs/genome_tool/gatk-4.6.2.0/gatk HaplotypeCaller -R {params.ref} '
        '--native-pair-hmm-threads {threads} --intervals chr16 '
        '-I {input} -D {params.dbsnp} '
        '-O {output} --tmp-dir /media/nvme8T/temp --ERC GVCF) 2>{log}'

rule HaplotypeCaller_chr17:
    input:
        '/media/nvme8T/WGS_pipeline_storage/{sample}.BQSR.soft_clip.hg38noalt.bam'
    output:
        temp('/media/nvme8T/WGS_pipeline/{sample}/{sample}.g.chr17.vcf.gz')
    params:
        ref = '/media/nvme1/reference/hg38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna',
        dbsnp = '/home/xs/genome_tool/reference-genome/hg38/dbsnp_146.hg38.vcf.gz'
    log:
        '/media/nvme8T/WGS_pipeline/{sample}/{sample}.chr17.HaplotypeCaller.log'
    resources: tmpdir='/media/nvme8T/temp'
    threads: 2
    shell:
        '(/home/xs/genome_tool/gatk-4.6.2.0/gatk HaplotypeCaller -R {params.ref} '
        '--native-pair-hmm-threads {threads} --intervals chr17 '
        '-I {input} -D {params.dbsnp} '
        '-O {output} --tmp-dir /media/nvme8T/temp --ERC GVCF) 2>{log}'

rule HaplotypeCaller_chr18:
    input:
        '/media/nvme8T/WGS_pipeline_storage/{sample}.BQSR.soft_clip.hg38noalt.bam'
    output:
        temp('/media/nvme8T/WGS_pipeline/{sample}/{sample}.g.chr18.vcf.gz')
    params:
        ref = '/media/nvme1/reference/hg38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna',
        dbsnp = '/home/xs/genome_tool/reference-genome/hg38/dbsnp_146.hg38.vcf.gz'
    log:
        '/media/nvme8T/WGS_pipeline/{sample}/{sample}.chr18.HaplotypeCaller.log'
    resources: tmpdir='/media/nvme8T/temp'
    threads: 2
    shell:
        '(/home/xs/genome_tool/gatk-4.6.2.0/gatk HaplotypeCaller -R {params.ref} '
        '--native-pair-hmm-threads {threads} --intervals chr18 '
        '-I {input} -D {params.dbsnp} '
        '-O {output} --tmp-dir /media/nvme8T/temp --ERC GVCF) 2>{log}'

rule HaplotypeCaller_chr19:
    input:
        '/media/nvme8T/WGS_pipeline_storage/{sample}.BQSR.soft_clip.hg38noalt.bam'
    output:
        temp('/media/nvme8T/WGS_pipeline/{sample}/{sample}.g.chr19.vcf.gz')
    params:
        ref = '/media/nvme1/reference/hg38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna',
        dbsnp = '/home/xs/genome_tool/reference-genome/hg38/dbsnp_146.hg38.vcf.gz'
    log:
        '/media/nvme8T/WGS_pipeline/{sample}/{sample}.chr19.HaplotypeCaller.log'
    resources: tmpdir='/media/nvme8T/temp'
    threads: 2
    shell:
        '(/home/xs/genome_tool/gatk-4.6.2.0/gatk HaplotypeCaller -R {params.ref} '
        '--native-pair-hmm-threads {threads} --intervals chr19 '
        '-I {input} -D {params.dbsnp} '
        '-O {output} --tmp-dir /media/nvme8T/temp --ERC GVCF) 2>{log}'

rule HaplotypeCaller_chr20:
    input:
        '/media/nvme8T/WGS_pipeline_storage/{sample}.BQSR.soft_clip.hg38noalt.bam'
    output:
        temp('/media/nvme8T/WGS_pipeline/{sample}/{sample}.g.chr20.vcf.gz')
    params:
        ref = '/media/nvme1/reference/hg38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna',
        dbsnp = '/home/xs/genome_tool/reference-genome/hg38/dbsnp_146.hg38.vcf.gz'
    log:
        '/media/nvme8T/WGS_pipeline/{sample}/{sample}.chr20.HaplotypeCaller.log'
    resources: tmpdir='/media/nvme8T/temp'
    threads: 2
    shell:
        '(/home/xs/genome_tool/gatk-4.6.2.0/gatk HaplotypeCaller -R {params.ref} '
        '--native-pair-hmm-threads {threads} --intervals chr20 '
        '-I {input} -D {params.dbsnp} '
        '-O {output} --tmp-dir /media/nvme8T/temp --ERC GVCF) 2>{log}'

rule HaplotypeCaller_chr21:
    input:
        '/media/nvme8T/WGS_pipeline_storage/{sample}.BQSR.soft_clip.hg38noalt.bam'
    output:
        temp('/media/nvme8T/WGS_pipeline/{sample}/{sample}.g.chr21.vcf.gz')
    params:
        ref = '/media/nvme1/reference/hg38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna',
        dbsnp = '/home/xs/genome_tool/reference-genome/hg38/dbsnp_146.hg38.vcf.gz'
    log:
        '/media/nvme8T/WGS_pipeline/{sample}/{sample}.chr21.HaplotypeCaller.log'
    resources: tmpdir='/media/nvme8T/temp'
    threads: 2
    shell:
        '(/home/xs/genome_tool/gatk-4.6.2.0/gatk HaplotypeCaller -R {params.ref} '
        '--native-pair-hmm-threads {threads} --intervals chr21 '
        '-I {input} -D {params.dbsnp} '
        '-O {output} --tmp-dir /media/nvme8T/temp --ERC GVCF) 2>{log}'

rule HaplotypeCaller_chr22:
    input:
        '/media/nvme8T/WGS_pipeline_storage/{sample}.BQSR.soft_clip.hg38noalt.bam'
    output:
        temp('/media/nvme8T/WGS_pipeline/{sample}/{sample}.g.chr22.vcf.gz')
    params:
        ref = '/media/nvme1/reference/hg38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna',
        dbsnp = '/home/xs/genome_tool/reference-genome/hg38/dbsnp_146.hg38.vcf.gz'
    log:
        '/media/nvme8T/WGS_pipeline/{sample}/{sample}.chr22.HaplotypeCaller.log'
    resources: tmpdir='/media/nvme8T/temp'
    threads: 2
    shell:
        '(/home/xs/genome_tool/gatk-4.6.2.0/gatk HaplotypeCaller -R {params.ref} '
        '--native-pair-hmm-threads {threads} --intervals chr22 '
        '-I {input} -D {params.dbsnp} '
        '-O {output} --tmp-dir /media/nvme8T/temp --ERC GVCF) 2>{log}'

rule HaplotypeCaller_chrM:
    input:
        '/media/nvme8T/WGS_pipeline_storage/{sample}.BQSR.soft_clip.hg38noalt.bam'
    output:
        temp('/media/nvme8T/WGS_pipeline/{sample}/{sample}.g.chrM.vcf.gz')
    params:
        ref = '/media/nvme1/reference/hg38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna',
        dbsnp = '/home/xs/genome_tool/reference-genome/hg38/dbsnp_146.hg38.vcf.gz'
    log:
        '/media/nvme8T/WGS_pipeline/{sample}/{sample}.chrM.HaplotypeCaller.log'
    resources: tmpdir='/media/nvme8T/temp'
    threads: 2
    shell:
        '(/home/xs/genome_tool/gatk-4.6.2.0/gatk HaplotypeCaller -R {params.ref} '
        '--native-pair-hmm-threads {threads} --intervals chrM '
        '-I {input} -D {params.dbsnp} '
        '-O {output} --tmp-dir /media/nvme8T/temp --ERC GVCF) 2>{log}'

rule HaplotypeCaller_merge:
    input:
        r1= '/media/nvme8T/WGS_pipeline/{sample}/{sample}.g.chrY.vcf.gz',
        r2= '/media/nvme8T/WGS_pipeline/{sample}/{sample}.g.chrX.vcf.gz',
        r3= '/media/nvme8T/WGS_pipeline/{sample}/{sample}.g.chr22.vcf.gz',
        r4= '/media/nvme8T/WGS_pipeline/{sample}/{sample}.g.chr21.vcf.gz',
        r5= '/media/nvme8T/WGS_pipeline/{sample}/{sample}.g.chr20.vcf.gz',
        r6= '/media/nvme8T/WGS_pipeline/{sample}/{sample}.g.chr19.vcf.gz',
        r7= '/media/nvme8T/WGS_pipeline/{sample}/{sample}.g.chr18.vcf.gz',
        r8= '/media/nvme8T/WGS_pipeline/{sample}/{sample}.g.chr17.vcf.gz',
        r9= '/media/nvme8T/WGS_pipeline/{sample}/{sample}.g.chr16.vcf.gz',
        r10= '/media/nvme8T/WGS_pipeline/{sample}/{sample}.g.chr15.vcf.gz',
        r11= '/media/nvme8T/WGS_pipeline/{sample}/{sample}.g.chr14.vcf.gz',
        r12= '/media/nvme8T/WGS_pipeline/{sample}/{sample}.g.chr13.vcf.gz',
        r13= '/media/nvme8T/WGS_pipeline/{sample}/{sample}.g.chr12.vcf.gz',
        r14= '/media/nvme8T/WGS_pipeline/{sample}/{sample}.g.chr11.vcf.gz',
        r15= '/media/nvme8T/WGS_pipeline/{sample}/{sample}.g.chr10.vcf.gz',
        r16= '/media/nvme8T/WGS_pipeline/{sample}/{sample}.g.chr9.vcf.gz',
        r17= '/media/nvme8T/WGS_pipeline/{sample}/{sample}.g.chr8.vcf.gz',
        r18= '/media/nvme8T/WGS_pipeline/{sample}/{sample}.g.chr7.vcf.gz',
        r19= '/media/nvme8T/WGS_pipeline/{sample}/{sample}.g.chr6.vcf.gz',
        r20= '/media/nvme8T/WGS_pipeline/{sample}/{sample}.g.chr5.vcf.gz',
        r21= '/media/nvme8T/WGS_pipeline/{sample}/{sample}.g.chr4.vcf.gz',
        r22= '/media/nvme8T/WGS_pipeline/{sample}/{sample}.g.chr3.vcf.gz',
        r23= '/media/nvme8T/WGS_pipeline/{sample}/{sample}.g.chr2.vcf.gz',
        r24= '/media/nvme8T/WGS_pipeline/{sample}/{sample}.g.chr1.vcf.gz',
        r25= '/media/nvme8T/WGS_pipeline/{sample}/{sample}.g.chrM.vcf.gz'
    output:
        '/media/nvme8T/WGS_pipeline_storage/{sample}.g.vcf.merge.gz'
    params:
        ref = '/media/nvme1/reference/hg38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna',
        dbsnp = '/home/xs/genome_tool/reference-genome/hg38/dbsnp_146.hg38.vcf.gz'
    log:
        '/media/nvme8T/WGS_pipeline/{sample}/{sample}.HaplotypeCaller.log'
    resources: tmpdir='/media/nvme8T/temp'
    shell:
        '(/home/xs/genome_tool/gatk-4.6.2.0/gatk MergeVcfs -I {input.r24} -I {input.r23} -I {input.r22} -I {input.r21} -I {input.r20} -I {input.r19} -I {input.r18} -I {input.r17} -I {input.r16} -I {input.r15} -I {input.r14} -I {input.r13} -I {input.r12} -I {input.r11} '
        '-I {input.r10} -I {input.r9} -I {input.r8} -I {input.r7} -I {input.r6} -I {input.r5} -I {input.r4} -I {input.r3} -I {input.r2} -I {input.r1} -I {input.r25} '
        '-O {output}) 2>{log}'

rule GenotypeGVCFs:
    input:
        '/media/nvme8T/WGS_pipeline_storage/{sample}.g.vcf.merge.gz'
    output:
        temp('/media/nvme8T/WGS_pipeline/{sample}/{sample}.vcf.gz')
    params:
        ref = '/media/nvme1/reference/hg38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna'
    log:
        '/media/nvme8T/WGS_pipeline/{sample}/{sample}.GenotypeGVCFs.log'
    resources: tmpdir='/media/nvme8T/temp'
    shell:
        '(/home/xs/genome_tool/gatk-4.6.2.0/gatk GenotypeGVCFs --tmp-dir /media/nvme8T/temp -R {params.ref} '
        '-V {input} -stand-call-conf 5 -O {output}) 2>{log}'

rule VariantRecalibrator_SNP:
    input:
        '/media/nvme8T/WGS_pipeline/{sample}/{sample}.vcf.gz'
    output:
        r1 = temp('/media/nvme8T/WGS_pipeline/{sample}/{sample}.VQSR.SNPS.recal'),
        r2 = temp('/media/nvme8T/WGS_pipeline/{sample}/{sample}.VQSR.snps.tranches')
    log:
        '/media/nvme8T/WGS_pipeline/{sample}/{sample}.VariantRecalibrator-SNP.log'
    params:
        ref = '/media/nvme1/reference/hg38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna',
        hapmap = ' /home/xs/genome_tool/reference-genome/hg38/hapmap_3.3.hg38.vcf.gz',
        omini = '/home/xs/genome_tool/reference-genome/hg38/1000G_omni2.5.hg38.vcf.gz',
        SNPS_1000G = '/home/xs/genome_tool/reference-genome/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz',
        dbsnp = '/home/xs/genome_tool/reference-genome/hg38/dbsnp_146.hg38.vcf.gz'
    resources: tmpdir='/media/nvme8T/temp'
    shell:
        '(/home/xs/genome_tool/gatk-4.6.2.0/gatk --java-options "-Xmx5G" VariantRecalibrator -R {params.ref} -V {input} '
        '-resource:hapmap,known=false,training=true,truth=true,prior=15.0 {params.hapmap} '
        '-resource:omini,known=false,training=true,truth=false,prior=12.0 {params.omini} '
        '-resource:1000G,known=false,training=true,truth=false,prior=10.0 {params.SNPS_1000G} '
        '-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {params.dbsnp} '
        '-an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum --max-gaussians 6 -mode SNP '
        '-tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0 '
        '--tranches-file {output.r2} --tmp-dir /media/nvme8T/temp '
        '-O {output.r1}) 2>{log}'
    
rule ApplyVQSR_SNP:
    input:
        r1 = '/media/nvme8T/WGS_pipeline/{sample}/{sample}.VQSR.SNPS.recal',
        r2 = '/media/nvme8T/WGS_pipeline/{sample}/{sample}.VQSR.snps.tranches',
        r3 = '/media/nvme8T/WGS_pipeline/{sample}/{sample}.vcf.gz'
    output:
        temp('/media/nvme8T/WGS_pipeline/{sample}/{sample}.VQSR.SNPs.vcf.gz')
    params:
        ref = '/media/nvme1/reference/hg38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna',
    log:
        '/media/nvme8T/WGS_pipeline/{sample}/{sample}.ApplyVQSR-SNP.log'
    resources: tmpdir='/media/nvme8T/temp'
    shell:
        '(/home/xs/genome_tool/gatk-4.6.2.0/gatk --java-options "-Xmx5G" ApplyVQSR --tmp-dir /media/nvme8T/temp -R {params.ref} -V {input.r3} '
        '-ts-filter-level 99.0 -mode SNP '
        '--tranches-file {input.r2} '
        '--recal-file {input.r1} '
        '-O {output}) 2>{log}'

rule VariantRecalibrator_indels:
    input:
        '/media/nvme8T/WGS_pipeline/{sample}/{sample}.VQSR.SNPs.vcf.gz'
    output:
        r1 = temp('/media/nvme8T/WGS_pipeline/{sample}/{sample}.VQSR.indels.recal'),
        r2 = temp('/media/nvme8T/WGS_pipeline/{sample}/{sample}.VQSR.indels.tranches')
    log:
        '/media/nvme8T/WGS_pipeline/{sample}/{sample}.VariantRecalibrator-indels.log'
    resources: tmpdir='/media/nvme8T/temp'
    params:
        ref = '/media/nvme1/reference/hg38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna',
        mills = '/home/xs/genome_tool/reference-genome/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz',
        dbsnp = '/home/xs/genome_tool/reference-genome/hg38/dbsnp_146.hg38.vcf.gz'
    shell:
        '(/home/xs/genome_tool/gatk-4.6.2.0/gatk --java-options "-Xmx5G" VariantRecalibrator --tmp-dir /media/nvme8T/temp -R {params.ref} -V {input} '
        '-resource:mills,known=true,training=true,truth=true,prior=12.0 {params.mills} '
        '-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {params.dbsnp} '
        '-an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum --max-gaussians 6 -mode INDEL '
        #'-an InbreedingCoeff '
        '-tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0 '
        '--tranches-file {output.r2} '
        '-O {output.r1}) 2>{log}'
    
rule ApplyVQSR_indels:
    input:
        r1 = '/media/nvme8T/WGS_pipeline/{sample}/{sample}.VQSR.indels.recal',
        r2 = '/media/nvme8T/WGS_pipeline/{sample}/{sample}.VQSR.indels.tranches',
        r3 = '/media/nvme8T/WGS_pipeline/{sample}/{sample}.VQSR.SNPs.vcf.gz'
    output:
        '/media/nvme8T/WGS_pipeline_storage/{sample}.VQSR.vcf.gz'
    params:
        ref = '/media/nvme1/reference/hg38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna'
    log:
        '/media/nvme8T/WGS_pipeline/{sample}/{sample}.ApplyVQSR-indels.log'
    resources: tmpdir='/media/nvme8T/temp'
    shell:
        '(/home/xs/genome_tool/gatk-4.6.2.0/gatk --java-options "-Xmx5G" ApplyVQSR --tmp-dir /media/nvme8T/temp -R {params.ref} -V {input.r3} '
        '-ts-filter-level 99.0 -mode INDEL '
        '--tranches-file {input.r2} '
        '--recal-file {input.r1} '
        '-O {output}) 2>{log}'

rule qualimap:
    input:
        r1 = '/media/nvme8T/WGS_pipeline_storage/{sample}.VQSR.vcf.gz',
        r2 = '/media/nvme8T/WGS_pipeline_storage/{sample}.BQSR.soft_clip.hg38noalt.bam'
    output:
        '/media/nvme8T/WGS_pipeline_storage/{sample}_qualimap_result/{sample}.pdf'
    params:
        dirname = '/media/nvme8T/WGS_pipeline_storage/{sample}_qualimap_result',
        outname = '{sample}.pdf'
    log:
        '/media/nvme8T/WGS_pipeline/{sample}/{sample}.qualimap.log'
    threads: 30
    shell:
        'mamba run -n base /home/xs/genome_tool/qualimap_v2.2.1/qualimap bamqc -bam {input.r2} '
        '-nr 3000 -outdir {params.dirname} -outfile {params.outname} -outformat PDF:HTML -nt {threads} '
        '--java-mem-size=60000M '  
