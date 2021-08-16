

configfile: 'config.yml'

ids = ['CRR034527']


rule all:
    input:
        expand('output/vartrix/{id}.mm', id=ids)


rule subset_vcf:
    input:
        vcf=config['dbsnp'],
        chrmap=config['chrmap'],
    output:
        'output/vcf/dbsnp_popfreq_chr6.vcf'
    shell:
        'bcftools view '
        '{input.vcf} '
        '-r NC_000006.12:1-58553888 | '
        'bcftools annotate --rename-chrs {input.chrmap} '
        '-o {output} '

rule vartrix:
    params:
        vartrix=config['vartrix'],
    input:
        vcf='output/vcf/dbsnp_popfreq_chr6.vcf',
        fasta=config['fasta'],
        bam=config['cellranger_base_path'] + '/{id}/outs/possorted_genome_bam.bam',
        barcodes=config['cellranger_base_path'] + '/{id}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz'
    output:
        'output/vartrix/{id}.mm'
    shell:
        '{params.vartrix} '
        '-v {input.vcf} '
        '-b {input.bam} '
        '-f {input.fasta} '
        '-c {input.barcodes} '
        '-o {output} '
