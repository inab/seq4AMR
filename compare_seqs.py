#!/usr/bin/python
import os
import Bio
import sys
import random
import hashlib
import argparse
import numpy as np
from Bio import SeqIO

def splitSequence(seq, length = 80):
  ''' Split a given sequence contained in one line into lines of size "length"
  '''
  return "\n".join([seq[i:i + length] for i in range(0, len(seq), length)])

if __name__ == "__main__":

  parser = argparse.ArgumentParser()

  parser.add_argument("-i", "--in", dest = "inFile", required = True, type = \
    str, help = "Input Sequence file")

  parser.add_argument("-r", "--ref", dest = "refFile", required = True, type = \
    str, help = "Input Reference Sequence file")

  parser.add_argument("-f", "--input_format", dest = "inFormat", type = str, \
    default = "fasta", help = "Set input alignment format")

  parser.add_argument("-o", "--out", dest = "outFile", default = None, type = \
    str, help = "Set output file")

  args = parser.parse_args()

  ## Check input parameters
  if not os.path.isfile(args.inFile):
    sys.exit(("ERROR: Check input sequences file '%s'") % (args.inFile))

  if not os.path.isfile(args.refFile):
    sys.exit(("ERROR: Check input sequences file '%s'") % (args.inFile))


  ## Read input sequences file and get some basic information from it e.g.
  ## sequences names, residues number, etc.
  reference = 0
  reference_seqs = {}
  for record in SeqIO.parse(args.refFile, args.inFormat):
    seq = str(record.seq)
    md5 = hashlib.md5(seq.encode('utf-8')).hexdigest()
    reference_seqs.setdefault(md5, {}).setdefault("seq", seq)
    reference_seqs.setdefault(md5, {}).setdefault("ids", []).append(record.id)
    reference += 1
  
  stats = {}
  input_seqs = {}
  for record in SeqIO.parse(args.inFile, args.inFormat):
    seq = str(record.seq)
    md5 = hashlib.md5(seq.encode('utf-8')).hexdigest()
    input_seqs.setdefault(md5, {}).setdefault("seq", seq)
    input_seqs.setdefault(md5, {}).setdefault("ids", []).append(record.id)

    if md5 in reference_seqs:
      stats.setdefault("hits", []).append(record.id)
    stats.setdefault("entries", []).append(record.id)
    
  ofile = open(args.outFile, "w") if args.outFile else sys.stdout

  ## Report true positive:
  positive = len(stats["hits"])
  all_entries = len(stats["entries"])
  true_positive_rate =  positive / reference
  accuracy = positive / all_entries
  
  print (("%-16s\t%d/%d\t%.4f") % ("True Positive Rate", positive, reference, \
    true_positive_rate))
  print (("%-16s\t%d/%d\t%.4f") % ("Accuracy", positive, all_entries, accuracy))
 
  ofile.close()
