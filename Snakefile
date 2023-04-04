

configfile: 'config.yml'

ids = config['ids']

output_path = 'output/' + config['version']

rule all:
    input:
        expand(output_path + '/cellsnp/{id}/cellSNP.tag.AD.mtx', id=ids),
        expand(        output_path + '/cellsnp-sces/{id}/cellsnp_sce_{id}.rds', id=ids),
        output_path + '/cellsnp-sce-list.rds',



# rule subset_vcf:
#     input:
#         vcf='snp-data/genome1K_sorted.vcf.gz',
#         chrmap=config['chrmap'],
#     output:
#         'output/vcf/genome1K_chr6.vcf'
#     shell:
#         'bcftools view '
#         '{input.vcf} '
#         '-r 6:1-29800000,6:33400000-58553888 | '
#         'bcftools annotate --rename-chrs {input.chrmap} '
#         '-o {output} '

rule subset_vcf2:
    input:
        vcf='snp-data/genome1K_sorted.vcf.gz',
        chrmap=config['chrmap'],
    output:
        output_path + '/vcf/genome1K_remapped.vcf'
    shell:
        'bcftools view '
        '{input.vcf} '
        '-r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 | '
        'bcftools annotate --rename-chrs {input.chrmap} '
        '-o {output} '

## Update March 21st 2023 to apply over all SNPs and not just 6p
rule cellsnp:
    params:
        cellsnp=config['cellsnp'],
    input:
        vcf=output_path + '/vcf/genome1K_remapped.vcf',
        bam=config['cellranger_base_path'] + '/{id}/outs/possorted_genome_bam.bam',
        barcodes=config['cellranger_base_path'] + '/{id}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz'
    output:
        output_path + '/cellsnp/{id}/cellSNP.tag.AD.mtx',
        output_path + '/cellsnp/{id}/cellSNP.tag.DP.mtx'
    shell:
        'mkdir -p output/cellsnp/{wildcards.id}/ && '
        '{params.cellsnp} '
        '-R {input.vcf} '
        '-s {input.bam} '
        '-b {input.barcodes} '
        '-O ' + output_path + '/cellsnp/{wildcards.id}/ '
        '--minMAF 0.1 --minCOUNT 20 '

rule cellsnp_to_sce:
    input:
        ad_mat = output_path + '/cellsnp/{id}/cellSNP.tag.AD.mtx',
        dp_mat = output_path + '/cellsnp/{id}/cellSNP.tag.DP.mtx'  ,
        sample = output_path + '/cellsnp/{id}/cellSNP.samples.tsv',
        snp_id = output_path + '/cellsnp/{id}/cellSNP.base.vcf',
    output:
        sce = output_path + '/cellsnp-sces/{id}/cellsnp_sce_{id}.rds'
    script:
        'cellsnp-to-sce.R'

rule cellsnp_sces_to_list:
    params:
        input_dir = output_path + '/cellsnp-sces/'
    input:
        expand(output_path + '/cellsnp-sces/{id}/cellsnp_sce_{id}.rds', id=ids)
    output:
        output_path + '/cellsnp-sce-list.rds'
    resources:
        mem_mb=10000
    shell:
        'Rscript sces-to-list.R '
        '--input_dir {params.input_dir} '
        '--output {output} '


rule to_expression_sces:
    params:
        dir = lambda wildcards: config['cellranger_base_path'] + f'/{wildcards.id}/outs/filtered_feature_bc_matrix/'
    input:
        config['cellranger_base_path'] + '/{id}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz'
    output:
        output_path + '/expression/sces/sce_{id}.rds'
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
    resources:
        mem_mb=10000
    shell:
        'Rscript concat-sce.R '
        '--input_dir output/expression/sces '
        '--output {output} '

# rule run_pileup_phase:
#     params:
#         samples = ','.join(ids),
#         bams = ','.join(expand(config['cellranger_base_path'] + '/{id}/outs/possorted_genome_bam.bam',id=ids)),
#         eagle = config['eagle'],
#         cellsnp = config['cellsnp'],
#         barcodes = ','.join(expand(config['cellranger_base_path'] + '/{id}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz',id=ids))
#     input:
#         bams = expand(config['cellranger_base_path'] + '/{id}/outs/possorted_genome_bam.bam',id=ids),
#         snpdata = config['snpdata'],
#         gmap = config['gmap']
#     output:
#         expand(output_path + 'pileup_phase/{id}_allele_counts.tsv.gz',
#         id=ids)
#     shell:
#         "mkdir -p output/phasing &&  "
#         "Rscript tools/pileup_and_phase.R  "
#         "--label test  "
#         "--samples {params.samples}  "
#         "--bams {params.bams} "
#         "--barcodes {params.barcodes} "
#         "--eagle {params.eagle} "
#         "--gmap {input.gmap} "
#         "--snpvcf {input.snpdata} "
#         "--cellsnp_lite {params.cellsnp} "
#         "--paneldir output/phasing  "
#         "--outdir output/pileup_phase "
#         "--ncores 1 "
