if (
    config["datatype"] == "SampleData[PairedEndSequencesWithQuality]"
    and config["jan-mode"] == False
):

    rule fastq_score:
        input:
            "results/{date}/out/joined-seqs.qza",
        output:
            filtering="results/{date}/out/demux-joined-filtered.qza",
            stats="results/{date}/out/demux-joined-filter-stats.qza",
        params:
            date=get_date(),
            min_quality=config["filtering"]["phred-score"],
            min_length_frac=config["filtering"]["min-length-frac"],
            max_ambig=config["filtering"]["max-ambiguity"],
        log:
            "logs/{date}/filtering/fastq-score.log",
        conda:
            "../envs/qiime-only-env.yaml"
        shell:
            "qiime quality-filter q-score "
            "--i-demux {input} "
            "--p-min-quality {params.min_quality} "
            "--p-min-length-fraction {params.min_length_frac} "
            "--p-max-ambiguous {params.max_ambig} "
            "--o-filtered-sequences {output.filtering} "
            "--o-filter-stats {output.stats} "
            "--verbose 2> {log}"


if (
    config["datatype"] == "SampleData[SequencesWithQuality]"
    and config["jan-mode"] == False
):

    rule fastq_score:
        input:
            "results/{date}/out/trimmed-seqs.qza",
        output:
            filtering="results/{date}/out/demux-joined-filtered.qza",
            stats="results/{date}/out/demux-joined-filter-stats.qza",
        params:
            date=get_date(),
            min_quality=config["filtering"]["phred-score"],
            min_length_frac=config["filtering"]["min-length-frac"],
            max_ambig=config["filtering"]["max-ambiguity"],
        log:
            "logs/{date}/filtering/fastq-score.log",
        conda:
            "../envs/qiime-only-env.yaml"
        shell:
            "qiime quality-filter q-score "
            "--i-demux {input} "
            "--p-min-quality {params.min_quality} "
            "--p-min-length-fraction {params.min_length_frac} "
            "--p-max-ambiguous {params.max_ambig} "
            "--o-filtered-sequences {output.filtering} "
            "--o-filter-stats {output.stats} "
            "--verbose 2> {log}"


if config["jan-mode"] == False:

    rule chimera_filtering:
        input:
            table="results/{date}/out/table-cluster.qza",
            seqs="results/{date}/out/seq-cluster.qza",
        output:
            direc=directory("results/{date}/out/uchime-dn-out"),
            table="results/{date}/out/table-nonchimeric-wo-borderline.qza",
            seqs="results/{date}/out/rep-seqs-nonchimeric-wo-borderline.qza",
        params:
            minh=config["filtering"]["chimera-minh"],
        log:
            "logs/{date}/filtering/chimera-filtering.log",
        conda:
            "../envs/qiime-only-env.yaml"
        shell:
            "qiime vsearch uchime-denovo "
            "--i-table {input.table} "
            "--i-sequences {input.seqs} "
            "--p-minh {params.minh} "
            "--output-dir {output.direc} \n"
            "qiime feature-table filter-features "
            "--i-table {input.table} "
            "--m-metadata-file {output.direc}/chimeras.qza "
            "--p-exclude-ids "
            "--o-filtered-table {output.table} \n"
            "qiime feature-table filter-seqs "
            "--i-data {input.seqs} "
            "--m-metadata-file {output.direc}/chimeras.qza "
            "--p-exclude-ids "
            "--o-filtered-data {output.seqs} "
            "--verbose 2> {log}"

    rule filter_seq_length:
        input:
            seq="results/{date}/out/rep-seqs-nonchimeric-wo-borderline.qza",  #results/{date}/out/seq-cluster.qza",
            table="results/{date}/out/table-nonchimeric-wo-borderline.qza",  #"results/{date}/out/table-cluster.qza"
        output:
            seq="results/{date}/out/seq-cluster-lengthfilter.qza",
            table="results/{date}/out/table-cluster-lengthfilter.qza",
        params:
            min_length=config["filtering"]["min-seq-length"],
        log:
            "logs/{date}/filtering/filter-seq-length.log",
        conda:
            "../envs/qiime-only-env.yaml"
        shell:
            "qiime feature-table filter-seqs "
            "--i-data {input.seq} "
            "--m-metadata-file {input.seq} "
            "--p-where 'length(sequence) > {params.min_length}' "
            "--o-filtered-data {output.seq} \n"
            "qiime feature-table filter-features "
            "--i-table {input.table} "
            "--m-metadata-file {output.seq} "
            "--o-filtered-table {output.table} "
            "--verbose 2> {log}"


if config["jan-mode"] == True:

    rule dada2:
        input:
            "results/{date}/out/trimmed-seqs.qza",
        output:
            table="results/{date}/out/table-cluster-lengthfilter.qza",
            seq="results/{date}/out/seq-cluster-lengthfilter.qza",
            stats="results/{date}/out/dada2-stats.qza",
        params:
            trunc_len_f=config["dada2"]["trunc-len-f"],
            trunc_len_r=config["dada2"]["trunc-len-r"],
            trim_left_f=config["dada2"]["trim-left-f"],
            trim_left_r=config["dada2"]["trim-left-r"],
            max_ee_f=config["dada2"]["max-ee-f"],
            max_ee_r=config["dada2"]["max-ee-r"],
            trunc_q=config["dada2"]["trunc-q"],
            min_overlap=config["dada2"]["min-overlap"],
            pooling_method=config["dada2"]["pooling-method"],
            chimera_method=config["dada2"]["chimera-method"],
            min_fold_parent_over_abundance=config["dada2"][
                "min-fold-parent-over-abundance"
            ],
            n_reads_learn=config["dada2"]["n-reads-learn"],
            threads=config["threads"],
        log:
            "logs/{date}/clustering/dada2.log",
        conda:
            "../envs/qiime-only-env.yaml"
        shell:
            "qiime dada2 denoise-paired "
            "--i-demultiplexed-seqs {input} "
            "--p-trunc-len-f {params.trunc_len_f} "
            "--p-trunc-len-r {params.trunc_len_r} "
            "--p-trim-left-f {params.trim_left_f} "
            "--p-trim-left-r  {params.trim_left_r} "
            "--p-max-ee-f {params.max_ee_f} "
            "--p-max-ee-r {params.max_ee_r} "
            "--p-trunc-q {params.trunc_q} "
            "--p-min-overlap {params.min_overlap} "
            "--p-pooling-method {params.pooling_method} "
            "--p-chimera-method {params.chimera_method} "
            "--p-min-fold-parent-over-abundance {params.min_fold_parent_over_abundance} "
            "--p-n-threads {params.threads} "
            "--p-n-reads-learn {params.n_reads_learn} "
            "--o-table {output.table} "
            "--o-representative-sequences {output.seq} "
            "--o-denoising-stats {output.stats} "
            "--verbose 2> {log}"


rule abundance_frequency:
    input:
        "results/{date}/visual/table-cluster-lengthfilter.qzv",
    output:
        abundance="results/{date}/out/abundance.txt",
        feature_table=report(
            directory("results/{date}/visual/table-cluster-lengthfilter/data"),
            caption="../report/feature-table.rst",
            category="4. Qualitycontrol",
            htmlindex="index.html",
        ),
    params:
        relative_abundance=config["filtering"]["relative-abundance-filter"],
    log:
        "logs/{date}/filtering/abundance-frequency.log",
    conda:
        "../envs/abundancefiltering.yaml"
    script:
        "../scripts/relative_abundance.py"


rule filter_frequency:
    input:
        table="results/{date}/out/table-cluster-lengthfilter.qza",  #"results/{date}/out/table-cluster.qza", 
        seqs="results/{date}/out/seq-cluster-lengthfilter.qza",  #"results/{date}/out/seq-cluster.qza", 
        abundance="results/{date}/out/abundance.txt",
    output:
        table="results/{date}/out/table-cluster-filtered.qza",  # "results/{date}/out/table-cluster-freq.qza"
        seqs="results/{date}/out/seq-cluster-filtered.qza",  # "results/{date}/out/seq-cluster-freq.qza"
    log:
        "logs/{date}/filtering/filter-frequency.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "value=$(<{input.abundance}) \n"
        "echo $value \n"
        "qiime feature-table filter-features "
        "--i-table {input.table} "
        "--p-min-frequency $value "
        "--o-filtered-table {output.table} \n"
        "qiime feature-table filter-seqs "
        "--i-data {input.seqs} "
        "--i-table {output.table} "
        "--p-no-exclude-ids "
        "--o-filtered-data {output.seqs}"


rule filter_human:
    input:
        seq="results/{date}/out/derepl-seq.qza",
        table="results/{date}/out/derepl-table.qza",
        ref_seq="resources/GRCh38_latest_genomic_upper.qza",
    output:
        seq="results/{date}/out/derep-seq-nonhum.qza",
        table="results/{date}/out/derep-table-nonhum.qza",
        human_hit="results/{date}/out/human.qza",
    params:
        threads=config["threads"],
        #perc-identity=config["filtering"]["perc-identity"],
        #perc-query-aligned=config["filtering"]["perc-query-aligned"],
    log:
        "logs/{date}/filtering/filter-human.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime quality-control exclude-seqs "
        "--i-query-sequences {input.seq} "
        "--i-reference-sequences {input.ref_seq} "
        "--p-threads {params.threads} "
        "--p-perc-identity 0.93 "
        "--p-perc-query-aligned 0.93 "
        "--o-sequence-hits {output.human_hit} "
        "--o-sequence-misses {output.seq} "
        "--verbose 2> {log} \n"
        "qiime feature-table filter-features "
        "--i-table {input.table} "
        "--m-metadata-file {output.seq} "
        "--o-filtered-table {output.table} "


rule taxa_collapse:
    input:
        table="results/{date}/out/table-cluster-filtered.qza",
        taxonomy="results/{date}/out/taxonomy.qza",
    output:
        "results/{date}/out/taxa_collapsed.qza",
    log:
        "logs/{date}/filtering/taxa-collapse.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime taxa collapse "
        "--i-table {input.table} "
        "--i-taxonomy {input.taxonomy} "
        "--p-level 6 "
        "--o-collapsed-table {output}"


rule filter_taxonomy:
    input:
        table="results/{date}/out/table-cluster-filtered.qza",
        seq="results/{date}/out/seq-cluster-filtered.qza",
        taxonomy="results/{date}/out/taxonomy.qza",
    output:
        table="results/{date}/out/table-taxa-filtered.qza",
        seq="results/{date}/out/seq-taxa-filtered.qza",
    log:
        "logs/{date}/filtering/filter-taxonomy.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime taxa filter-table "
        "--i-table {input.table} "
        "--i-taxonomy {input.taxonomy} "
        "--p-exclude mitochondria,chloroplast "
        "--o-filtered-table {output.table} \n"
        "qiime taxa filter-seqs "
        "--i-sequences {input.seq} "
        "--i-taxonomy {input.taxonomy} "
        "--p-exclude mitochondria,chloroplast "
        "--o-filtered-sequences {output.seq} "
