proportional_htseq
==================

This script, designed to work under Linux environments, will help you to run htseq-count taking the multihits into account proportionally. 
Each read will be weighted according to the number of mapped locations. For example, a read mapped to 5 different positions will add 0.2 to the counts of each feature.

Dependencies:
- samtools 0.1.18 or above (http://samtools.sourceforge.net/)
- htseq 0.5.3p3 or above (http://www-huber.embl.de/users/anders/HTSeq/doc/count.html)


Usage: python htseq.py [options] \<mandatory>

    Options:
        -h, --help:
                 show this help message and exit
        -m, --mode:
                 mode to handle reads overlapping more than one feature(choices: union, intersection-strict, intersection-nonempty; default: union)
        -s, --stranded:
                 whether the data is from a strand-specific assay. Specify 'yes', 'no', or 'reverse' (default: yes). 'reverse' means 'yes' with reversed strand interpretation
        -a, --minaqual:
                 skip all reads with alignment quality lower than the given minimum value (default: 0)
        -t, --type:
                 feature type (3rd column in GFF file) to be used, all features of other type are ignored (default, suitable for Ensembl GTF files: exon)
        -i, --idattr:
                 GFF attribute to be used as feature ID (default, suitable for Ensembl GTF files: gene_id)
        -o, --samout:
                 write out all SAM alignment records into an output SAM file called SAMOUT, annotating each line with its feature assignment (as an optional field with tag 'XF')
        -n, --sort:
                 sort the bam file by name (necessary for paired-end reads
        -d, --destiny:
                 output directory (default, the directory of execution)
        -c, --clean_up:
                 Do not remove the intermediate files generated (default: remove intermediate files)
        -p, --threads:
                 Number of threads to run (default: 1)
    Mandatory:
        -b, --bam:
                 bam/sam file to read
        -g, --gtf:
                 gtf file

    07/03/2013. Pedro Furió.
