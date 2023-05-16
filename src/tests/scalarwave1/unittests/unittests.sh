TOP=../../../..

NMESH=$TOP/exe/nmesh
FDIFF=$TOP/utilities/floatdiff

$NMESH sw7FV-DG.par > sw7FV-DG.out
$FDIFF -at 1e-13 -rt 1e-12 sw7FV-DG_correct sw7FV-DG
