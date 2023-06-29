TOP=../../../..

rm -rf *_correct *_previous

NMESH=$TOP/exe/nmesh
mpirun -np 2 $NMESH misc_unit1.par	> misc_unit1.out

mv -v misc_unit1			misc_unit1_correct
