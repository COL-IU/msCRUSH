""" Parse clusters formed by SpectraCluster software.
    
    Args:
        inFname: file name of the result of SpectraCluster, for e.g.
            'Adult_Liver.mgf.spectracluster'.
        ouFname: output file containing clusters of spectra.
    Returns:
        A file containing clusters formed by SpectraCluster software. 
        The result file comes with a header 'ID Titles', then follows each line
        a cluster. For example

        ID Titles 
        3e68f113-5fbf-48ec-a75b-de0902641e76    Adult_Liver_Gel_Elite_83_f05.2507.2507.2.dta;Adult_Liver_Gel_Elite_83_f02.4256.4256.2.dta;Adult_Liver_Gel_Elite_83_f02.5351.5351.2.dta;

        You can find example output file in ../../data/Adult_Liver_cluster_spectracluster
"""


import sys
import time


inFname = sys.argv[1]
ouFname = sys.argv[2]
input_title = sys.argv[3]
charge = int(sys.argv[4])

CLUSTER = '=Cluster='
SPECTRUM = 'SPEC'
ID = 'id'
TITLE_DELIM = 'title='
CONSENSUS_MZ = 'consensus_mz'
CONSENSUS_INTEN = 'consensus_intens'
PRECURSORMZ = "av_precursor_mz"

mzs = []
intens = []
precursor_mz = ''
spectrum_idx = 0
with open(ouFname, 'w') as ouf:
    with open(inFname) as inf:
        for line in inf:
          line = line.rstrip()
          if PRECURSORMZ in line:
            precursor_mz = line.split('=')[1]

          if CONSENSUS_MZ in line:
            mzs = line.split('=')[1].split(',')

          if CONSENSUS_INTEN in line:
            intens = line.split('=')[1].split(',')

            spectrum_idx = spectrum_idx + 1
            # Write spectrum to file.
            ouf.write("BEGIN IONS\n")
            ouf.write("PEPMASS={0}\n".format(precursor_mz))
            ouf.write("CHARGE={0}+\n".format(charge))
            #ouf.write("TITLE=PXD001197.{0}.{0}.{1}.dta\n".format(spectrum_idx, charge))
            ouf.write("TITLE={2}.{0}.{0}.{1}.dta\n".format(spectrum_idx, charge, input_title))
            for mz, inten in zip(mzs, intens):
              ouf.write("{0}\t{1}\n".format(mz, inten))
            ouf.write("END IONS\n")
            ouf.write("\n")

print ("congrats! spectra extraction from SpectraCluster results is finished!")

