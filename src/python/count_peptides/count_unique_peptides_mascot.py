import argparse
import numpy as np
import pandas as pd
import sys, time

FILESCAN = 'FileScan'
SCORE = 'Score'
SEQUENCE = 'Sequence'
DECOY = 'Decoy'
REFERENCE = 'Reference'
REVERSED = 'REVERSED'

FLAGS=None

def main():
  mascot_file = FLAGS.mascot_file
  print 'mascot file: ', mascot_file
  title2peptide = {}
  peptides = set()
  decoy_peptides = set()

# Read in mascot result.
  df = pd.read_csv(mascot_file, sep='\t')
#print list(df)

  min_score_for_decoy = 9999
  decoy_cnt = 0
  for i in range(len(df)):

    if df[FILESCAN][i] is not np.nan:
      # Check for decoy.
      peptide = df[SEQUENCE][i]
      peptide = peptide.split('.')[1]
      title = df[FILESCAN][i]

      # Remove modifications.
      if FLAGS.remove_ptm:
        peptide = peptide.replace('*', '')

      # Transform I->L.
      if FLAGS.enable_i2l:
        peptide = peptide.replace('I', 'L')

      # Convert I amino acid to L amino acid.
      # peptide = peptide.replace('I', 'L')

      # Change to another way to decide decoy or not, so that compatible with Windows platform result.
      if REVERSED in df[REFERENCE][i]:
        decoy_cnt = decoy_cnt + 1
        if peptide not in decoy_peptides:
          decoy_peptides.add(peptide)
      else:

        if title in title2peptide:
          print 'found duplicates', title
          continue

        title2peptide[title] = peptide

        if peptide not in peptides:
          peptides.add(peptide)
    else:
      pass
      #print 'found exception at line: ',i

  print '#decoy psm count: ', decoy_cnt
  print '#decoy peptides: ', len(decoy_peptides)
  print '#unique peptides:', len(peptides)
  print 'title2peptide dict size:', len(title2peptide)

  if FLAGS.output_peptides:
    output_file = FLAGS.mascot_file + '.mascot_peptides'

    if FLAGS.enable_i2l:
      output_file += '.i2l'

    if FLAGS.remove_ptm:
      output_file += '.rm-ptm'

    print 'Writing peptides to file', output_file

    with open(output_file, 'w') as of:
      of.write('\n'.join(peptides))
      of.write('\n')
    print 'write complete.'
  else:
    print 'print a random unique peptide: ', peptides.pop()

  if FLAGS.output_psm:
    output_file = FLAGS.mascot_file + '.psm'

    if FLAGS.enable_i2l:
      output_file += '.i2l'
    if FLAGS.remove_ptm:
      output_file += '.rm-ptm'

    print 'Writing peptides to file', output_file

    with open(output_file, 'w') as of:
      of.write('Title\tPeptide\n')
      for title in title2peptide:
        of.write('{0}\t{1}\n'.format(title, title2peptide[title]))
    print 'write complete.'

if __name__ == '__main__':
  parser = argparse.ArgumentParser()

  parser.add_argument('--mascot_file', type=str,
                      default=None,
                      help='mascot file to parse.')

  parser.add_argument('--output_peptides', type=bool,
                      default=False,
                      help='Whether to ouput unique peptides.')

  parser.add_argument('--output_psm', type=bool,
                      default=False,
                      help='Whether to ouput PSM.')

  parser.add_argument('--enable_i2l', type=bool,
                      default=False,
                      help='Whether to enable peptide seq I->L.')

  parser.add_argument('--remove_ptm', type=bool,
                      default=False,
                      help='Whether to remove PTM in peptide seq.')

  FLAGS, unparsed = parser.parse_known_args()
  main()

