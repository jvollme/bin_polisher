## po2pocp.py
Purify pre-generated bins (e.g. using Maxbin) based on z-score differences in contig coverage using multiple sequencing datasets 

````
optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show version info and quit

Required input arguments:
  -if INPUT_FASTA_LIST [INPUT_FASTA_LIST ...], --input_fasta INPUT_FASTA_LIST [INPUT_FASTA_LIST ...]
                        Input fasta file(s) of a (single) bin (may be compressed)
  -ic INPUT_COVERAGE_LIST [INPUT_COVERAGE_LIST ...], --input_coverage INPUT_COVERAGE_LIST [INPUT_COVERAGE_LIST ...]
                        seperate abundance files for each dataset, listing the
                        respective abundance/coverage of each contig in a
                        seperate line (Format: <contig-id>\TAB<coverage>), May
                        include contigs not present in bin-fasta (such contigs
                        wil be ignored)

Filtering options:
  -pr {high,low,both,none}, --pre-remove {high,low,both,none}
                        remove extreme values based on upper (99%) or lower
                        (1%) percentile, or both. default = none
  -mi MAX_ITERATIONS, --max_iterations MAX_ITERATIONS
                        maximum number of iterations for recalculating and
                        comparing coverage average, standard deviation and
                        z-score values and removing difference-outliers.
                        Default = 50
  -uzc UPPER_ZSCORE_CUTOFF, --upper_zscore_cutoff UPPER_ZSCORE_CUTOFF
                        Start iteratively removing contigs with z-score
                        differences above cutoff starting at the specified
                        value (cutoff will be decreased by 1 when no z-scores
                        above cutoff are encountered, until the lower zscore
                        cutoff is reached). Default = 4
  -lzc LOWER_ZSCORE_CUTOFF, --lower_zscore_cutoff LOWER_ZSCORE_CUTOFF
                        Stop iteratively reducing z-score difference cutoff if
                        it falls below this value. Default = 2

Output options:
  --intermediate        output results of all intermediate iterations. Default
                        = only output results of last iteration
  --out_bad             create seperate output files for rejected reads also.
                        default = False
  -op OUT_PREFIX, --out_prefix OUT_PREFIX
                        prefix for output file(s). default = "bin_polisher"
````

example usage: ```bin_polisher.py -if <bin1.fasta.gz> -ic <bin1_sample1.coverage> [<bin1_sample2.coverage> ...]```

This tool was originally created for lab-internal use and not with user friendliness or adaptiveness in mind. It is published here, to enable reproducability of any of our research projects that may have been using this script (e.g. [Vollmers et al, 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5372823/)).

More flexible versions may be uploaded whenever I have time for it. For any problems or questions about the usage, please create an Issue on this repository.


