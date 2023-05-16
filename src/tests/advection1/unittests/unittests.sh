TOP=../../../..

NMESH=$TOP/exe/nmesh
FDIFF=$TOP/utilities/floatdiff

$NMESH advDG-FD.par > advDG-FD.out
$FDIFF -at 1e-13 -rt 1e-12 advDG-FD_correct advDG-FD
