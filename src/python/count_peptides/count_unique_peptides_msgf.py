""" Count the unique peptides in an MSGF MS identification results.
    Args:
        inFname: input file name of MSGF+ identification results.        
        charge: charge of MS.

    Return:
        
    Example:
        python count_unique_peptides_msgf.py XXX
         
"""
import argparse
import numpy as np
import pandas as pd
import re
import sys
import time

FLAGS = None

CHARGE = 'Charge'
DECOY = 'XXX_'
# MSGFSCORE = 'MSGFScore'
MSGFSCORE = 'DeNovoScore'
PEPTIDE = 'Peptide'
PROTEIN = 'Protein'
PEPQVALUE = 'PepQValue'
QVALUE = 'QValue'
SPECID = 'SpecID'
TITLE = 'Title'

delim = '\t'

def main():
  inFname = FLAGS.msgf_file
  charge = FLAGS.charge


  print 'file: ', inFname
  print 'charge: ', charge
  if FLAGS.enable_i2l:
    print 'enable i->l'

  # Remove PTM
  if FLAGS.remove_ptm:
    print 'remove ptm'

  # Convert PTM to star (*)
  if FLAGS.ptm_2_star:
    print 'remove ptm->star'

  df = pd.read_csv(inFname, sep=delim)

  #print list(df.columns.values)
  #print FLAGS.use_qvalue

  if FLAGS.use_qvalue:
    standard = QVALUE
    print 'Use QVALUE for FDR.'
  else:
    standard = PEPQVALUE
    print 'Use PEPQVALUE for FDR.'

  # Needed use the regex more than once, so better to be compiled.
  if FLAGS.ptm_2_star:
    pattern = re.compile(r"\+?\d+\.?\d+")

  title2peptide = {}
  peptides = set()
  peptide2MSGFScore = {}
  peptide2MSGFScoreAvg = {}
  for i in range(0, len(df)):
      protein = df[PROTEIN][i]
      # Get rid of reversed peptides.
      if DECOY in protein:
          # print(protein)
          continue

      if charge != int(df[CHARGE][i]):
          continue

      if float(df[standard][i]) > 0.01:
          continue

      # Check TITLE in the column, because later we use different msgf converter
      if TITLE in df.columns:
        title = df[TITLE][i]
      elif SPECID in df.columns:
        title = df[SPECID][i]
      else:
        print 'TITLE not correct.'
        sys.exit(-1)

      peptide = df[PEPTIDE][i]
      # Remove N- and C- terminal.
      if '.' == peptide[1] and '.' == peptide[-2]:
        peptide = peptide[2:][:-2]

      # Keep only need the first PSM.
      if title in title2peptide:
          #print 'found duplicates', title
          #sys.exit(-1)
          continue

      # Transform I->L.
      if FLAGS.enable_i2l:
        peptide = peptide.replace('I', 'L')

      # Remove PTM
      if FLAGS.remove_ptm:
        peptide = peptide.translate(None, '0123456789+.')

      # Convert PTM to star (*)
      if FLAGS.ptm_2_star:
        peptide = pattern.sub("*", peptide)

      title2peptide[title] = peptide

      msgfScore = float(df[MSGFSCORE][i])
      if peptide not in peptide2MSGFScore:
        peptide2MSGFScore[peptide] = msgfScore
      elif msgfScore > peptide2MSGFScore[peptide]:
        peptide2MSGFScore[peptide] = msgfScore

      if peptide not in peptide2MSGFScoreAvg:
        peptide2MSGFScoreAvg[peptide] = [msgfScore]
      else:
        peptide2MSGFScoreAvg[peptide].append(msgfScore)

      
      #print title, peptide,  msgfScore
      #time.sleep(1)

      if peptide not in peptides:
        peptides.add(peptide)

  print 'MSGF Peptide-level Max Score Mean: ', np.mean(peptide2MSGFScore.values())
  print 'MSGF Peptide-level Max Score std: ', np.std(peptide2MSGFScore.values())

  print 'MSGF Peptide-level Max Score min: ', min(peptide2MSGFScore.items(), key=lambda x: x[1])
  print 'MSGF Peptide-level Max Score max: ', max(peptide2MSGFScore.items(), key=lambda x: x[1])

  scores = [np.mean(peptide2MSGFScoreAvg[peptide]) for peptide in peptide2MSGFScoreAvg]
  print 'MSGF Peptide-level Avg Score Mean: ', np.mean(scores)
  print 'MSGF Peptide-level Avg Score std: ', np.std(scores)


  print 'identified #unique peptides: ', len(peptides)
  print 'identified #PSM(1 PSM/spectrum): ', len(title2peptide)

  if FLAGS.output_peptides:
    output_file = FLAGS.msgf_file + '.msgf_peptides'
    
    if FLAGS.enable_i2l:
      output_file += '.i2l'
    if FLAGS.remove_ptm:
      output_file += '.rm-ptm'
    if FLAGS.ptm_2_star:
      output_file += '.ptm2star'

    print 'Writing peptides to file', output_file
    with open(output_file, 'w') as of:
      of.write('\n'.join(peptides))
      of.write('\n')

  if FLAGS.output_psm:
    output_file = FLAGS.msgf_file + '.psm'

    if FLAGS.enable_i2l:
      output_file += '.i2l'
    if FLAGS.remove_ptm:
      output_file += '.rm-ptm'
    if FLAGS.ptm_2_star:
      output_file += '.ptm2star'

    print 'Writing PSM to file', output_file

    with open(output_file, 'w') as of:
      of.write('Title\tPeptide\n')
      for title in title2peptide:
        of.write('{0}\t{1}\n'.format(title, title2peptide[title]))

  print 'program finished successfully!\n\n'


if __name__ == '__main__':
  parser = argparse.ArgumentParser()

  parser.add_argument('--msgf_file', type=str,
                      default=None,
                      help='msgf file to analyze.')

  parser.add_argument('--charge', type=int,
                      default=-1,
                      help='charge of MS in msgf file.')

  parser.add_argument('--output_peptides', type=bool,
                      default=False,
                      help='Whether to ouput unique peptides.')

  parser.add_argument('--output_psm', type=bool,
                      default=False,
                      help='Whether to ouput PSM.')

  parser.add_argument('--use_qvalue', type=bool,
                      default=False,
                      help='Whether to QValue to do fdr.')

  parser.add_argument('--enable_i2l', type=bool,
                      default=False,
                      help='Whether to enable peptide seq I->L.')

  parser.add_argument('--remove_ptm', type=bool,
                      default=False,
                      help='Whether to remove PTM in peptide seq.')

  parser.add_argument('--ptm_2_star', type=bool,
                      default=False,
                      help='Whether to convert PTM to star(*).')

  FLAGS, unparsed = parser.parse_known_args()
  main()

