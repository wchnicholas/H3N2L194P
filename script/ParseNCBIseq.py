#!/usr/bin/python
import os
import re
import sys
import glob
from Bio import SeqIO
from collections import Counter, defaultdict

def isInt(astring):
    """ Is the given string an integer? """
    try: int(astring)
    except ValueError: return 0
    else: return 1

def main():
  P194Lpos = 209
  alnfile  = 'Fasta/Bris07_fromNCBI.aln'
  records  = [record for record in SeqIO.parse(alnfile,"fasta")]
  for record in records:
    ID = str(record.id)
    resi194 = str(record.seq)[P194Lpos]
    print ID+"\t"+resi194

if __name__ == "__main__":
  main()
