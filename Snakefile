

configfile: 'config.yml'

ids = ['CRR034499',
    'CRR034500',
    'CRR034501',
    'CRR034503',
    'CRR034504',
    'CRR034506',
    'CRR034507',
    'CRR034509',
    # 'CRR034510',
    'CRR034513',
    'CRR034516',
    'CRR034517',
    'CRR034518',
    'CRR034520',
    'CRR034524',
    'CRR034527']


rule all:
    input:
        expand('output/cellsnp/{id}/cellSNP.tag.AD.mtx', id=ids)


rule subset_vcf:
    input:
        vcf='snp-data/genome1K_sorted.vcf.gz',
        chrmap=config['chrmap'],
    output:
        'output/vcf/genome1K_chr6.vcf'
    shell:
        'bcftools view '
        '{input.vcf} '
        '-r 6:1-58553888 | '
        'bcftools annotate --rename-chrs {input.chrmap} '
        '-o {output} '

rule cellsnp:
    params:
        cellsnp=config['cellsnp'],
    input:
        vcf='output/vcf/genome1K_chr6.vcf',
        bam=config['cellranger_base_path'] + '/{id}/outs/possorted_genome_bam.bam',
        barcodes=config['cellranger_base_path'] + '/{id}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz'
    output:
        'output/cellsnp/{id}/cellSNP.tag.AD.mtx',
        'output/cellsnp/{id}/cellSNP.tag.DP.mtx'
    shell:
        'mkdir -p output/cellsnp/{wildcards.id}/ && '
        '{params.cellsnp} '
        '-R {input.vcf} '
        '-s {input.bam} '
        '-b {input.barcodes} '
        '-O output/cellsnp/{wildcards.id}/ '
        '--minMAF 0.1 --minCOUNT 20 '


# ../../software/cellsnp-lite/cellsnp-lite -s /home/campbell/share/datasets/peng-2019-cellranger/CRR034527/outs/possorted_genome_bam.bam -b /home/campbell/share/datasets/peng-2019-cellranger/CRR034527/outs/filtered_feature_bc_matrix/barcodes.tsv.gz -O test/ -p 4 --minMAF 0.1 --minCOUNT 20 -R output/vcf/genome1K_chr6.vcf