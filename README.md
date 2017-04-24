## H3N2 egg-adaptive substitution L194P
This analysis is adapted from [McWhite et al. 2016](https://academic.oup.com/ve/article-lookup/doi/10.1093/ve/vew026)

### Input File
* Fasta/pdmH1N1\_All.fa: 2009 pandemic H1N1 (swine flu) HA sequences downloaded from [GISAID](http://platform.gisaid.org)
* Fasta/HumanH3N2\_All.fa: Human H3N2 HA sequences downloaded from [GISAID](http://platform.gisaid.org)
  * Since there is a limit on the number of sequences being downloaded at once on GISAID, sequences for this project was first downloaded separately based on continent. Then sequences from different continent were combined to a single Fasta file.
* Fasta/Bris07\_fromNCBI.fa: 11 Bris07 sequences from the NCBI protein database (https://www.ncbi.nlm.nih.gov/protein/) were obtained by searching “A/Brisbane/10/2007, hemagglutinin”. 

### Protocol
1. Multiple sequence alignment using [MAFFT version 7.157b](http://mafft.cbrc.jp/alignment/software/)
* mafft --auto Fasta/pdmH1N1\_All.fa > Fasta/pdmH1N1\_All.aln
* mafft --auto Fasta/HumanH3N2\_All.fa > Fasta/HumanH3N2\_All.aln
* mafft --auto Fasta/Bris07\_fromNCBI.fa > Fasta/Bris07\_fromNCBI.aln
