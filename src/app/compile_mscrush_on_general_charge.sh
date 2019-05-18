res=mscrush_on_general_charge
echo 'compiling msCRUSH clustering program for general charge...'
g++ -std=c++11 -O3 -D FAST_PARSE ../utility/io.cc ../class/distance.cc  ../class/core.cc greedyCluster_omp_general_charge.cc -o $res -fopenmp 
#echo 'exe:' $res
