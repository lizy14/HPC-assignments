methods=( "atomic" "critical_extra_var" "critical" "reduction" "reduction_dynamic" "serial" )
for i in "${methods[@]}"
do
    echo $i
    cp $i.h header.h
    g++ source.cpp -o exe -std=c++11 -fopenmp -DNON_VERBOSE 
    ./exe
done
