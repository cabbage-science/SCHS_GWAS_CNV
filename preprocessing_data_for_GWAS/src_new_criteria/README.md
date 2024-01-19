Source files for processing of raw PennCNV output, after taking into account criteria used by Ryan Collin's 741k papeer
- Filtering of samples in the PHENOTYPE file, remove those that "LRR_SD > 0.3, BAF_drift > 0.01, WF > |0.02| and CNV > 100" - did not fulfil PennCNV QC
- removing of rows from PennCNV output in blacklisted regions
- Creation of new matrix according to 2 new criteria of CNVs (remember need multiple scripts because there is too much data)
- Splitting of new matrix into individual chromosomes
