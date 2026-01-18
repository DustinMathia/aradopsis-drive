rule aggregate_hybrids:
    input: 
        # Points to samples/{sample}.csv
        f"{config['samples_dir']}/{{sample}}.csv"
    output: 
        # Points to results/{sample}.aggregate.csv
        f"{config['data_dir']}/{{sample}}.aggregate.csv"
    script:
        "../scripts/count_hybrids.py" # Note: Script path is relative to the Snakefile location

rule make_plots:
    input:
        f"{config['data_dir']}/{{sample}}.aggregate.csv"
    output:
        f"{config['results_dir']}/{{sample}}.h_ratio.pdf"
    script:
        "../scripts/plot_h_ratio.R"