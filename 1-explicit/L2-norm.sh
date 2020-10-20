#!/bin/bash
# Shell script to generate pair-wise L² norm at *t*=200

echo -e "Pairwise\n"

echo "dx1,dx2,l2"
./mmsp2norm run-d/refine-2/nuc81d2.80000.dat    run-d/refine-1/nuc81d1.20000.dat
./mmsp2norm run-d/refine-3/nuc81d3.320000.dat   run-d/refine-2/nuc81d2.80000.dat
./mmsp2norm run-d/refine-4/nuc81d4.1280000.dat  run-d/refine-3/nuc81d3.320000.dat
#./mmsp2norm run-d/refine-5/nuc81d5.5120000.dat  run-d/refine-4/nuc81d4.1280000.dat
#./mmsp2norm run-d/refine-6/nuc81d6.20480000.dat run-d/refine-5/nuc81d5.5120000.dat
#./mmsp2norm run-d/refine-7/nuc81d7.81920000.dat run-d/refine-6/nuc81d6.20480000.dat

#echo -e "\nMost refined\n"
#EXEMPLAR=run-d/run-d/refine-6/nuc81d6.20480000.dat

#echo "dx1,dx2,l2"
#./mmsp2norm ${EXEMPLAR} run-d/refine-1/nuc81d1.20000.dat
#./mmsp2norm ${EXEMPLAR} run-d/refine-2/nuc81d2.80000.dat
#./mmsp2norm ${EXEMPLAR} run-d/refine-3/nuc81d3.320000.dat
#./mmsp2norm ${EXEMPLAR} run-d/refine-4/nuc81d4.1280000.dat
#./mmsp2norm ${EXEMPLAR} run-d/refine-5/nuc81d5.5120000.dat
#./mmsp2norm ${EXEMPLAR} run-d/refine-6/nuc81d6.20480000.dat
