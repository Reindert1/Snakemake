#!/usr/bin/env python3

configfile: "config/config.yaml"

rule all:
    """ final rule """
    input: "result/heatmap.jpg"


rule make_heatmap:
    """ makes heatmap from data"""
    input: config["data"]

    output: "result/heatmap.jpg"

    run:
        from snakemake.utils import R
        R("""
        d <- as.matrix(read.csv('{input}', header=FALSE, sep=",")[-1,-1])
        rownames(d) <- read.csv('{input}', header=FALSE, sep=",")[-1,1]
        colnames(d) <- read.csv('{input}', header=FALSE, sep=",")[1,-1] 
        
        jpeg("{output}")
        heatmap(d) 
        
        dev.off()
        """)









