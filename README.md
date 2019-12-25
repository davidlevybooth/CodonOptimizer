# CodonOptimizer
Codon optimization script built on DNAChisel with RSCU calculation from the CAI module

# Installation
Simply download codonusage.py or clone repository

# Usage
python codonusage.py <options>
  
-i  --input_fasta  Path to the gene to be optimized (fasta file) (Required)

-e  --highexp  Path to input genes to calculate codon usage (RSCU) (Optional)

-t  --taxid  NCBI Taxonomic ID # (e.g., 101510 for Rhodococcus jostii RHA1 (either -e or it are Required)

-c  --genetic_code  Genetic code # currently only code 11 (Bacteria, Archaea) is supported (Optional)

-r  --report  Path to report output file. Must have a .zip extention

-o  --output  Path to write output codon optimized fasta file (Required)

-m  --method  Method of codon optimization (see below) Default = match_codon_usage

  
# Method Decription
Codon-optimize a coding sequence using a user-selected method. This pseudo-specification is actually a function which returns an instance of another specification class depending on the selected "method":
    - For method="use_best_codon", every codon will be replaced by the "best"
      (i.e. most frequent) synonymous codon in the target organism. This is
      equivalent to Codon Adaptation Index (CAI) optimization.
    - For method="match_codon_usage", the final sequence's codon usage will
      match as much as possible the codon usage profile of the target species
      (this method is used throughout the literature, see for instance Hale
      and Thomson 1998).
    - For method="harmonize_rca", Each codon will be replaced by a synonymous
      codon whose usage in the target organism matches the usage of the
      original codon in its host organism (as per Claassens 2017).
      
# Description
A simple script to wrap the codon optimization tools in the DNA Chisel package. Use codonusage.py to optimize the codon usage by either matching the codon usage of an input nucelotide fasta file, or from a known taxonomic id number (See https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi). An example sequence set for highly expressed genes from Rhodococcus jostii RHA1 is provided (RHA1_high_expression.fasta). The output can take two forms: 1) only the optimized sequence or 2) the optimized sequence and a brief DNA Chisel report. 

