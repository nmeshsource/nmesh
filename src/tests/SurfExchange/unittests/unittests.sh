TOP=../../../..

NMESH=$TOP/exe/nmesh
FDIFF=$TOP/utilities/floatdiff

/bin/cp SurfExchange_unit1.par SurfExchange_unit1_p1.par
/bin/cp SurfExchange_unit1.par SurfExchange_unit1_p2.par
/bin/cp SurfExchange_unit1.par SurfExchange_unit1_p3.par

mpirun -np 1 $NMESH SurfExchange_unit1_p1.par > SurfExchange_unit1.out1
mpirun -np 2 $NMESH SurfExchange_unit1_p2.par > SurfExchange_unit1.out2
mpirun -np 3 $NMESH SurfExchange_unit1_p2.par > SurfExchange_unit1.out3

$FDIFF -at 1e-13 -rt 1e-12 SurfExchange_unit1_p1_corr SurfExchange_unit1_p1
$FDIFF -at 1e-13 -rt 1e-12 SurfExchange_unit1_p2_corr SurfExchange_unit1_p2
$FDIFF -at 1e-13 -rt 1e-12 SurfExchange_unit1_p3_corr SurfExchange_unit1_p3
