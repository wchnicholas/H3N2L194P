#!/usr/bin/python
import os
import re
import sys
import glob
from scipy import stats
from Bio import SeqIO
from collections import Counter, defaultdict

def isInt(astring):
    """ Is the given string an integer? """
    try: int(astring)
    except ValueError: return 0
    else: return 1

def ExtractEggPSGNumber(PSG):
  if '/' not in PSG and '+' not in PSG and ',' not in PSG: 
    PSG_number = PSG.replace('_','')
    if ':' in PSG_number: PSG_number = PSG_number.rsplit(':')[1]
    if 'E' in PSG_number and isInt(PSG_number[1::]):
      return PSG_number.replace('E','')
    elif 'AM' in PSG_number and isInt(PSG_number[2::]):
      return PSG_number.replace('AM','')

def ParseSeqs(records,seq_outfile,pos):
  outfile  = open(seq_outfile,'w')
  Egg_dict = defaultdict(list)
  Ori_dict = defaultdict(list)
  excludepattern = re.compile ("UNKNOWN_1_RHMK|TMK1_MDCK|AMNIOTIC_1_PRHMK_2|M1_RII1,C5|R1_C|R1_S|RII1_C|RII1_S|RIIX_C|RX_C|MDCK_1_RHMK|NC|_MK1")
  unpassagedpattern = re.compile("LUNG|P0|OR_|ORIGINAL|CLINICAL|DIRECT")
  eggpattern = re.compile("AM[1-9]|E[1-7]|AMNIOTIC|EGG|EX|AM_[1-9]")
  cellpattern = re.compile("S[1-9]|SX|SIAT|MDCK|C[1-9]|CX|C_[1-9]|M[1-9]|MX|X[1-9]|^X_$")
  siatpattern = re.compile("^S[1-9]_$|SIAT2_SIAT1|SIAT3_SIAT1")
  monkeypattern=re.compile("TMK|RMK|RHMK|RII|PMK|R[1-9]|RX")
  siatexcludepattern=re.compile("SIAT|SX|S[1-9]")
  CountSeq = 0
  Egg_PSGNumber_dict = defaultdict(list)
  for record in records:
    CountSeq += 1
    header = str(record.id)
    if header.count('|')!=4: continue
    ID   = header.rsplit('|')[0]
    PSG  = header.rsplit('|')[1]
    year = header.rsplit('|')[-1][1:5]
    seq  = str(record.seq)
    assert(isInt(year))
    if 'X' in seq: continue
    if excludepattern.search(PSG): continue
    elif eggpattern.search(PSG): 
      Egg_dict[year].append(seq)
      PSG_Number = ExtractEggPSGNumber(PSG)
      if PSG_Number: Egg_PSGNumber_dict['All'].append(PSG_Number)
      if seq[pos]=='P':
        if PSG_Number: Egg_PSGNumber_dict['P'].append(PSG_Number)
        outfile.write('>Egg_P194_'+str(CountSeq)+"\n"+seq.replace('-','')+"\n")
      else:
        outfile.write('>Egg'+str(CountSeq)+"\n"+seq.replace('-','')+"\n")
    elif unpassagedpattern.search(PSG): 
      Ori_dict[year].append(seq)
      outfile.write('>Ori'+str(CountSeq)+"\n"+seq.replace('-','')+"\n")
  outfile.close()
  return Egg_dict, Ori_dict, Egg_PSGNumber_dict
    
def ExtractPosOfInterest(seq_dict,pos,PSG,outfile):
  for year in sorted(seq_dict.keys(),key=lambda x:int(x)):
    aa_dict = defaultdict(int)
    for seq in seq_dict[year]:
      aa = seq[pos]
      aa_dict[aa]+=1
    for aa in aa_dict.keys():
      outfile.write("\t".join(map(str,[year, PSG, aa, aa_dict[aa]]))+"\n")

def ComplilePSGNumber(Egg_PSGNumber_dict, PSG_outfile):
  outfile  = open(PSG_outfile, 'w') 
  P_dict   = Counter(Egg_PSGNumber_dict['P'])
  All_dict = Counter(Egg_PSGNumber_dict['All'])
  PSGnumbers = sorted(list(set(P_dict.keys()+All_dict.keys())))
  outfile.write("\t".join(['Passage Number', 'Pro', 'All', 'Proportion of Pro','SE'])+"\n")
  for PSGnumber in PSGnumbers:
    P_count   = P_dict[PSGnumber]
    All_count = All_dict[PSGnumber]
    P_frac    = float(P_count)/float(All_count)
    P_SE      = (P_frac*(1-P_frac)/All_count)**0.5
    outfile.write("\t".join(map(str,[PSGnumber, P_count, All_count, P_frac, P_SE]))+"\n")
  outfile.close()

def wrapper(alnfilename, P194Lpos, outfilename, seqoutfilename, PSGoutfilename):
  alnfile  = alnfilename
  outfile  = open(outfilename,'w')
  seq_outfile = seqoutfilename
  PSG_outfile = PSGoutfilename
  records  = [record for record in SeqIO.parse(alnfile,"fasta")]
  Egg_dict, Ori_dict, Egg_PSGNumber_dict = ParseSeqs(records,seq_outfile,P194Lpos)
  outfile.write("\t".join(['Year', 'Passage', 'AA', 'Count'])+"\n")
  ExtractPosOfInterest(Egg_dict,P194Lpos,'Egg',outfile)
  ExtractPosOfInterest(Ori_dict,P194Lpos,'Ori',outfile)
  outfile.close()
  if 'H3N2' in alnfile: ComplilePSGNumber(Egg_PSGNumber_dict, PSG_outfile)

def main():
  wrapper('Fasta/HumanH3N2_All.aln',225,'result/HumanH3N2_Pos194YearVsPSG.tsv',
          'result/HumanH3N2_EggOri.fa','result/HumanH3N2_PSG.tsv')
  wrapper('Fasta/pdmH1N1_All.aln',216,'result/pdmH1N1_Pos194YearVsPSG.tsv',
          'result/pdmH1N1_EggOri.fa','result/pdmH1N1_PSG.tsv')

if __name__ == "__main__":
  main()
