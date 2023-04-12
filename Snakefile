shell.executable("/bin/bash")

# Snakemake pipeline for extracting and clustering HsdS genes from xylella genomes.

#---- Fixed Variables ----#
configfile: "config/pipeline_config.yaml"

#genomes_dir = config["genomes_dir"]
genomes_dir = "genomes"
IDS, = glob_wildcards(genomes_dir + "/{sample}.fna")
# HsdQueries, = glob_wildcards("queries/blast/Riv5.{HsdQuery}.fna")
HSDX = [ "hsdR", "hsdM", "hsdS" ]
RMS = ["XfaI", "XfaII", "XfaIII", "XfaIV"]

#---- Pipeline ----#


rule all:
    input:
        expand("01-blastdbs/{sample}.fna.ndb", sample = IDS),
        expand("02-blastn/{Hsd}/{sample}.{Hsd}.{component}.blast.csv", sample = IDS, Hsd = RMS, component = HSDX),
        expand("03-blastdbcmd/{RMS}/{sample}.{RMS}.{HSDX}.fna", sample = IDS, RMS = RMS, HSDX = HSDX),
        expand("05-CD_HIT/inputs/{RMS}.hsdS.hits.and.refs.fna", RMS = RMS),
        expand("05-CD_HIT/output/{RMS}/{RMS}.clusters.clstr", RMS = RMS),
        expand("05-CD_HIT/output/{RMS}/multi-seq/0", RMS = RMS)
        # expand("04-hits_summary/concat.{RMS}.{HSDX}.hits.csv", sample = IDS, RMS = RMS, HSDX = HSDX)
# #        expand("03-blastdbcmd/", sample = IDS)


# rule makeblastdb:
    # input: "genomes/{sample}.fna"
# #    #expand(["{dir}/{sample}.{ext}"], dir = genomes_dir, ext = config["extension"], sample = IDS)    # each strain in genomes/*.fna
    # output:
        # nhdfile = "01-blastdbs/{sample}.nhd"
    # params:
        # dbname = "01-blastdbs/{sample}"
    # threads: 
        # workflow.cores
    # shell:
        # "makeblastdb -in {sample} -out {params.dbname} -dbtype nucl -hash_index -parse_seqids"
rule blast_makedatabase_nucleotide:
    input:
        fasta="genomes/{sample}.fna"
    output:
        multiext("01-blastdbs/{sample}.fna",
            ".ndb",
            ".nhr",
            ".nin",
            ".not",
            ".nsq",
            ".ntf",
            ".nto"
        )
    log:
        "logs/{sample}.log"
    params:
        "-input_type fasta -blastdb_version 5 -parse_seqids -hash_index"
    wrapper:
        "v1.25.0/bio/blast/makeblastdb"

#        query = expand("queries/blast/Riv5.{iRMS}.{iHSDX}.fna", iRMS = RMS, iHSDX = HSDX),

rule blast_nucleotide:
    input:
        query=("queries/blast/Riv5.{iRMS}.{iHSDX}.fna"),
        blastdb=multiext("01-blastdbs/{sample}.fna",
            ".ndb",
            ".nhr",
            ".nin",
            ".not",
            ".nsq",
            ".ntf",
            ".nto"
        )
    output:
        "02-blastn/{iRMS}/{sample}.{iRMS}.{iHSDX}.blast.csv"
    log:
        "logs/{sample}.{iRMS}.{iHSDX}.blast.log"
    threads:
        workflow.cores
    params:
        # Usable options and specifiers for the different output formats are listed here:
        # https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/blast/blastn.html.
        format="10 std slen",
        extra= expand("-evalue {param} -max_target_seqs 1", param = config["evalue"]) 
    wrapper:
        "v1.25.0/bio/blast/blastn"

rule blastdbcmd:
    input: 
        "02-blastn/{RMS}/{sample}.{RMS}.{HSDX}.blast.csv",
#    "{rules.blast_nucleotide.output}"
    output:
        "03-blastdbcmd/{RMS}/{sample}.{RMS}.{HSDX}.fna"
    params:
        blastdb="01-blastdbs/{sample}.fna"
    threads:
        workflow.cores
    run: 
        shell("[ -s {input} ] && blastdbcmd -db {params.blastdb} -entry $(cat {input} | sort -t',' -k12 -n -r | head -n1 | cut -f2 -d',') -out {output} || touch {output}"),
        shell('sed -i "s/>/>{wildcards.sample} /g" {output}') 
    
    

# rule blast_hits_summary:
    # input: 
        # expand("03-blastdbcmd/{RMS}/{sample}.{RMS}.{HSDX}.fna", sample = IDS, RMS = RMS, HSDX = HSDX)
    # output:
        # "04-hits_summary/concat.{RMS}.{HSDX}.hits.csv"
    # shell: 
        # "echo {wildcards.RMS},{wildcards.sample},$(head -n1 {input}) >> {output}" # make a script snippet?

rule concatenate_hsdS_hits:
    input:
        expand("03-blastdbcmd/{RM}/{sample}.{RM}.hsdS.fna", sample = IDS, RM= RMS)
    output:
        "05-CD_HIT/inputs/{RM}.hsdS.hits.and.refs.fna"
    params:
        catfiles = "03-blastdbcmd/{RM}/*.{RM}.hsdS.fna",
        refalleles = "queries/reference_hsdS/{RM}_type_alleles_04142022.fas"
    shell: 
        "cat {params.catfiles} {params.refalleles} > {output}"

rule cd_hit:
    input: 
        "05-CD_HIT/inputs/{RMS}.hsdS.hits.and.refs.fna"
    output: 
        fasta = "05-CD_HIT/output/{RMS}/{RMS}.clusters",
        clusters = "05-CD_HIT/output/{RMS}/{RMS}.clusters.clstr"
    # params:
        # outprefix = "05-CD_HIT/output/{RMS}/{RMS}.clusters"
    threads: 
        workflow.cores
    shell:
        "cd-hit-est -i {input} -o {output.fasta} -c 0.9 -aL 0.9 -d 0 -g 1 -sc 1"
        
rule make_multi_seq:
    input: 
        "05-CD_HIT/output/{RMS}/{RMS}.clusters.clstr" # cd-hit output
    output: 
        "05-CD_HIT/output/{RMS}/multi-seq/0"
    params:
        fastas = "05-CD_HIT/inputs/{RMS}.hsdS.hits.and.refs.fna",
        outpath = "05-CD_HIT/output/{RMS}/multi-seq/"
    run:
        shell("make_multi_seq.pl {params.fastas} {input} {params.outpath}")
        

        
        
        
        
        
# unused rules
# rule blastn:
    # input: 
        # db = "{rule.makeblastdb.output.nhdfile}",
        # query = "queries/blast/Riv5.{HsdQuery}.fna"
    # output:
        # "02-blastn/{sample}.{HsdQuery}.hits.csv"
    # params:
        # dbname = "01-blastdb/{sample}",
        # evalue = expand("{param}", param = config["evalue"])
    # threads:
        # workflow.cores
    # shell:
        # "blastn -db {params.dbname} -outfmt '10 std slen' -evalue {params.evalue} -out {output} -query {input.query} -max_target_seqs 1 -num_threads {threads}"


