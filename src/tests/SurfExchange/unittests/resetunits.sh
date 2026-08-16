TOP=../../../..

rm -rf *_correct *_previous

NMESH=$TOP/exe/nmesh
mpirun -np 2 $NMESH SurfExchange_unit1.par	> SurfExchange_unit1.out

mv -v SurfExchange_unit1			SurfExchange_unit1_correct
