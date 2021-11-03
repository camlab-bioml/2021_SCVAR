

configfile: 'config.yml'

ids = config['ids']


rule all:
    input:
        expand('output/cellsnp/{id}/cellSNP.tag.AD.mtx', id=ids),
        expand('output/expression/sces/sce_{id}.rds', id=ids),
        'output/expression/sce.rds'


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


rule to_expression_sces:
    params:
        dir = lambda wildcards: config['cellranger_base_path'] + f'/{wildcards.id}/outs/filtered_feature_bc_matrix/'
    input:
        config['cellranger_base_path'] + '/{id}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz'
    output:
        'output/expression/sces/sce_{id}.rds'
    shell:
        'Rscript to-expression-sce.R '
        '--input_dir {params.dir} '
        '--id {wildcards.id} '
        '--output {output} '

rule to_expression_sce:
    input:
        expand('output/expression/sces/sce_{id}.rds', id=ids),
    output:
        'output/expression/sce.rds'
    shell:
        'Rscript concat-sce.R '
        '--input_dir output/expression/sces '
        '--output {output} '