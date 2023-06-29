TOP=../../../..

NMESH=$TOP/exe/nmesh
FDIFF=$TOP/utilities/floatdiff

mpirun -np 2 $NMESH misc_unit1.par > misc_unit1.out
$FDIFF -at 1e-13 -rt 1e-12 misc_unit1_correct misc_unit1
