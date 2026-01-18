rule aggregate_hybrids:
    input: 
        # Points to samples/{sample}.csv
        f"{config['samples_dir']}/{{sample}}.csv"
    output: 
        # Points to results/{sample}.aggregate.csv
        f"{config['results_dir']}/{{sample}}.aggregate.csv"
    script:
        "../scripts/count_hybrids.py" # Note: Script path is relative to the Snakefile location