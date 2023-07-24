if config["bowtie"] == True:
    #rule create_bowtie_db:
    #    input:
    #        "resources/GRCh38_latest_genomic_upper.fna"
    #    output:
    #        directory("resources/bowtie_host_DB/")
    #    log:
    #        "logs/bowtie_db.log"
    #    conda:
    #        "../envs/python.yaml"
    #    shell:
    #        "bowtie2-build {input} {output}"

    rule map_sequences:
        input:
            db="resources/GRCh38_noalt_as",
            read1="data/{date}/{names}_L001_R1_001.fastq.gz",
            read2="data/{date}/{names}_L001_R2_001.fastq.gz",
        output:
            "results/{date}/out/bowtie/{names}_mapped_and_unmapped.sam"
        log:
            "logs/{date}/bowtie/{names}_mapping.log"
        conda:
            "../envs/python.yaml"
        shell:
            "bowtie2 -p 8 -x {input.db}/GRCh38_noalt_as "
            "-1 {input.read1} "
            "-2 {input.read2} "
            "--un-conc-gz "
            "SAMPLE_host_removed > {output} "
            "2> {log} "

    rule sam_to_bam:
        input:
            "results/{date}/out/bowtie/{names}_mapped_and_unmapped.sam"
        output:
            "results/{date}/out/bowtie/{names}_mapped_and_unmapped.bam"
        log:
            "logs/{date}/bowtie/{names}_sam_to_bam.log"
        conda:
            "../envs/python.yaml"
        shell:
            "samtools view -bS {input} > {output} "
            "2> {log} " 

    rule filter_unmapped:
        input:
            "results/{date}/out/bowtie/{names}_mapped_and_unmapped.bam"
        output:
            "results/{date}/out/bowtie/{names}_bothReadsUnmapped.bam"
        log:
            "logs/{date}/bowtie/{names}_filter_unmapped.log"
        conda:
            "../envs/python.yaml"
        shell:
            "samtools view -b -f 12 -F 256 "
            "{input} > {output} "
            "2> {log} "

    rule split_paired:
        input:
            #expand(
            "results/{date}/out/bowtie/{names}_bothReadsUnmapped.bam",
            #    names=get_reads_for_kraken(),
            #),
        output:
            #expand(
            read1="results/{date}/bowtie/data/{names}_L001_R1_001.fastq.gz",
            read2="results/{date}/bowtie/data/{names}_L001_R2_001.fastq.gz",
            sorte="results/{date}/bowtie/{names}_bothReadsUnmapped_sorted.bam",
             #   names=get_reads_for_kraken(),
            #),
            #dir="results/{date}/bowtie/data/",
        #params:
        # sorted="results/{date}/bowtie/{names}_bothReadsUnmapped_sorted.bam
        log:
            "logs/{{date}}/bowtie/{names}_split_paired.log"
        conda:
            "../envs/python.yaml"
        shell:
            "samtools sort -n -m 5G -@ 2 {input} -o {output.sorte} 2> {log} \n"
            "samtools fastq -@ 8 {output.sorte} "
            "-1 {output.read1} "
            "-2 {output.read2} "
            "-0 /dev/null -s /dev/null -n "
            "2> {log} "

    rule get_human_reads:
        input:
            "results/{date}/out/bowtie/{names}_mapped_and_unmapped.bam"
        output:
            "results/{date}/out/bowtie/readcount/{names}_aligned_reads_count.txt"
        log:
            "logs/{date}/bowtie/{names}_humanreads.log"
        conda:
            "../envs/python.yaml"
        shell:
            "samtools view -F 4 -c {input} > {output} "
            "2> {log} "


    rule complete_human_count:
        input:
            expand(
                "results/{{date}}/out/bowtie/readcount/{names}_aligned_reads_count.txt",
                names=get_reads_for_kraken(),
            )
        output:
            report(
                "results/{date}/visual/sample_frequencys_difference.csv",
                caption="../report/hum_filter_difference.rst",
                category="4. Qualitycontrol",
            ),
        log:
            "logs/{date}/bowtie/humanreads_complete.log"
        conda:
            "../envs/python.yaml"
        script:
            "../scripts/bowtie_frequency.py"

