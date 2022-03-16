rule dereplication:
    input:
        "results/{date}/out/demux-joined-filtered.qza",
    output:
        table="results/{date}/out/derepl-table.qza",
        seqs="results/{date}/out/derepl-seq.qza",
    shell:
        "qiime vsearch dereplicate-sequences "
        "--i-sequences {input} "
        "--o-dereplicated-table {output.table} "
        "--o-dereplicated-sequences {output.seqs}"


rule de_novo_clustering:
    input:
        table="results/{date}/out/derepl-table.qza",
        seqs="results/{date}/out/derepl-seq.qza",
    output:
        table="results/{date}/out/table-cluster.qza",
        seqs="results/{date}/out/seq-cluster.qza",
    params:
        perc_identity=0.99,
    shell:
        "qiime vsearch cluster-features-de-novo "
        "--i-table {input.table} "
        "--i-sequences {input.seqs} "
        "--p-perc-identity {params.perc_identity} "
        "--p-threads 10 "
        "--o-clustered-table {output.table} "
        "--o-clustered-sequences {output.seqs}"


rule classification:
    input:
        query="results/{date}/out/seq-cluster-filtered.qza",
        reference_reads=config["database"]["sequences"],
        reference_taxonomy=config["database"]["taxonomy"],
    output:
        "results/{date}/out/taxonomy.qza",
    params:
        perc_identity=0.97,
    shell:
        "qiime feature-classifier classify-consensus-vsearch "
        "--i-query {input.query} "
        "--i-reference-reads {input.reference_reads} "
        "--i-reference-taxonomy {input.reference_taxonomy} "
        "--p-maxaccepts 1 "
        "--p-maxrejects 1 "
        "--p-perc-identity {params.perc_identity} "
        "--p-threads 10 "
        "--o-classification {output} "
        "--verbose"


rule phylogenetic_tree:
    input:
        "results/{date}/out/seq-taxa-filtered.qza",  #"results/{date}/out/seq-cluster-filtered.qza"
    output:
        alignment="results/{date}/out/aligned-rep-seqs.qza",
        masked_alignment="results/{date}/out/masked-aligned-rep-seqs.qza",
        tree="results/{date}/visual/unrooted-tree.qza",
        rooted_tree="results/{date}/visual/rooted-tree.qza",
    shell:
        "qiime phylogeny align-to-tree-mafft-fasttree "
        "--i-sequences {input} "
        "--o-alignment {output.alignment} "
        "--o-masked-alignment {output.masked_alignment} "
        "--o-tree {output.tree} "
        "--o-rooted-tree {output.rooted_tree}"


rule core_metrics:
    input:
        phylogeny="results/{date}/visual/rooted-tree.qza",
        table="results/{date}/table-cluster-filtered.qza",
        metadata="config/pep/sample.tsv",
    output:
        "results/{date}/core-metrics-results",
    params:
        depth=100,
    shell:
        "qiime diversity core-metrics-phylogenetic "
        "--i-phylogeny {input.phylogeny} "
        "--i-table {input.table} "
        "--p-sampling-depth {params.depth} "
        "--m-metadata-file {input.metadata} "
        "--output-dir {output}"
