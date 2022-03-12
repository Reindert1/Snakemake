datadir= '/commons/Themas/Thema11/Dataprocessing/WC04/data/'

rule all:
    """ final rule """
    input: 'result/histogram.jpg'


rule make_histogram:
    """ rule that creates histogram from gene expression counts"""
    input:
        datadir + 'out.csv'
    output:
         'result/histogram.jpg'
    run:
        from snakemake.utils import R
        R("""
        data <-read.csv(file = "{input}", header=FALSE, sep=";")
        jpeg("{output}")
        hist(x=data$V1,
             main="expression values of gene CCND3 Cyclin D3",
             ylab="count",
             xlab="gene CCND3 Cyclin D3 expression")
        dev.off()
        """)
