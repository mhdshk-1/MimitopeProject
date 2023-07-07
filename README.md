# Summary
These scripts constitute a pipeline simulating MHC-II antigen processing. Given the input of a (viral) proteome, we use (Netchop)[https://services.healthtech.dtu.dk/services/NetChop-3.1/] to predict protein cleavage, carrying the resulting peptides through for binding prediction with (NetMHC2pan)[https://github.com/GfellerLab/MixMHC2pred]. Strongly-binding putative epitopes resulting from these steps are then compared against another proteome to identify potentially cross-reactive epitopes. This pipeline is presently a work in progress, and is not yet fully functional. Ideally we need to reach sufficient sensitivity such that we are able to detect the cross-reactive viral epitopes and proteins, and their associated HLA risk alleles, described in docs/worksheet.docx (corresponding to the following autoimmune diseases: Type 1 Diabetes, Narcolepsy Type 1, and Multiple Sclerosis). The below descriptions are summaries; see the scripts themselves for more detailed information.

# Netchop Scripts
1. netchop-1.py chunks the starting *.fasta proteome into small enough files for Netchop webserver submission.
2. netchop-2.py processes the outputs and combines them into a single *.csv file.

# MixMHC2pred Scripts
1. mixmhc2pred-1.py coerces the *.csv file into a *.tsv file for MixMHC2pred submission.
2. mixmhc2pred-2.py submits this file via command line along with the requested HLA alleles for binding prediction, assuming the program is accessable on system PATH.
3. mixmhc2pred-3.py processes the outputs for all alleles and combines them into a single *.csv file, categorising strong and weak binders.

# Cross-Reactivity Scripts
sequence-alignment-netchop-input.py take the *.csv file of strong binders and compare them against another proteome, to identify potentially cross-reactive epitopes. The match sequence length and identity level are adjustable parameters.