res=generate_consensus_spectrum_for_mscrush
echo 'compiling consensus spectra generation program...'
g++ -std=c++11 -O3 -D FAST_PARSE ../utility/io.cc ../class/core.cc  ../class/distance.cc generate_consensus_spectrum_for_mscrush.cc -o $res
#echo 'exe:' $res
