TOP=../../../..

NMESH=$TOP/exe/nmesh
FDIFF=$TOP/utilities/floatdiff

mpirun -np 2 $NMESH SurfExchange_unit1.par > SurfExchange_unit1.out
$FDIFF -at 1e-13 -rt 1e-12 SurfExchange_unit1_correct SurfExchange_unit1
