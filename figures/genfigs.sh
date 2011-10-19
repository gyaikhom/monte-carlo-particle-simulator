#!/bin/bash

make;
for x in `ls *.mps`; do
    f=`basename $x .mps`;
    ./mpstoeps.sh $x;
    epstopdf $f.eps;
done;
make clean;
