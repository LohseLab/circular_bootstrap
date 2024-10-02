***Usage***
```
mamba install -c conda-forge -c bioconda numpy pandas matplotlib
./circular_bootstrap.py input.tsv
```
Assumes the input formatted as follows with windows sorted in sequential order on each chromosome:

Chromosome  window_label  feature_1  ..  feature_n

... 


The presence of a header is assumed but the headers themselves are not specific.
