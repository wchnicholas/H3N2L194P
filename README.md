This README describes the scripts used for the sequence analysis in:   
[A structural explanation for the low effectiveness of the seasonal influenza H3N2 vaccine](http://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1006682)

## H3N2 egg-adaptive substitution L194P
This analysis is adapted from [McWhite et al. 2016](https://academic.oup.com/ve/article-lookup/doi/10.1093/ve/vew026)

### Input File
* Fasta/pdmH1N1\_All.fa: 2009 pandemic H1N1 (swine flu) HA sequences downloaded from [GISAID](http://platform.gisaid.org)
* Fasta/HumanH3N2\_All.fa: Human H3N2 HA sequences downloaded from [GISAID](http://platform.gisaid.org)
  * Since there is a limit on the number of sequences being downloaded at once on GISAID, sequences for this project was first downloaded separately based on continent. Then sequences from different continent were combined to a single Fasta file.
* Fasta/Bris07\_fromNCBI.fa: 11 Bris07 sequences from the [NCBI protein database](https://www.ncbi.nlm.nih.gov/protein/) were obtained by searching **"A/Brisbane/10/2007", hemagglutinin**.
* Fasta/HK14\_fromGISAID.fa: 8 HK14 sequences from [GISAID](http://platform.gisaid.org) were obtained by **"A/Hong Kong/4801/2014"**.
* Fasta/Sing16\_fromGISAID.fa: 4 Sing16 sequences from [GISAID](http://platform.gisaid.org) were obtained by searching **"A/Singapore/INFIMH-16-0019/2016"**.

### Protocol
#### Analyzing HA sequences from egg-passaged isolates downloaded from GISAID
1. Multiple sequence alignment (MSA) using [MAFFT version 7.157b](http://mafft.cbrc.jp/alignment/software/)
* mafft --auto Fasta/pdmH1N1\_All.fa > Fasta/pdmH1N1\_All.aln
* mafft --auto Fasta/HumanH3N2\_All.fa > Fasta/HumanH3N2\_All.aln

2. Parse MSA files to extract information on egg-passaged isolates
* python script/ParseGISAIDaln.py: 
  * Input files:
    * Fasta/pdmH1N1\_All.aln
    * Fasta/HumanH3N2\_All.aln
  * Output files:
    * result/HumanH3N2\_Pos194YearVsPSG.tsv
    * result/HumanH3N2\_EggOri.fa
    * result/HumanH3N2\_PSG.tsv
    * result/pdmH1N1\_Pos194YearVsPSG.tsv
    * result/pdmH1N1\_EggOri.fa
    * result/pdmH1N1\_PSG.tsv

3. Plot the frequency of different amino acids observed at residue 194 in different year
* Rscript script/Plot\_YearVsPSG.R
  * Input files:
    * result/H3N2\_Pos194YearVsPSG.tsv
    * result/pdmH1N1\_Pos194YearVsPSG.tsv
  * Output files: 
    * graph/H3N2\_YearVsAA\_resi194.png
    * graph/pdmH1N1\_YearVsAA\_resi194.png

4. Plot the frequency of L194P against the number of passage in eggs
* Rscript script/Plot\_ProVsPSG.R
  * Input file: 
    * result/HumanH3N2\_PSG.tsv
  * Output file:
    * graph/HumanH3N2\_ProVsPSG.png

#### Analyzing HA sequences from A/Brisbane/10/2007 (H3N2) deposited in NCBI
1. Multiple sequence alignment (MSA) using [MAFFT version 7.157b](http://mafft.cbrc.jp/alignment/software/) 
* mafft --auto Fasta/Bris07\_fromNCBI.fa > Fasta/Bris07\_fromNCBI.aln

2. Parse MSA files to extract amino-acid identity on residue 194
* python script/ParseNCBIseq.py
  * Input file: 
    * Fasta/Bris07\_fromNCBI.aln
  * Standard Output
