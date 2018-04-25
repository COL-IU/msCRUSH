res=splash_on_general_charge
echo 'compiling SPLASH clustering on general charge...'
g++ -std=c++11 -O3 ../utility/io.cc ../class/distance.cc  ../class/core.cc greedyCluster_omp_general_charge.cc -o $res -fopenmp 
echo 'exe:' $res
