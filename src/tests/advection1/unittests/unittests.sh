TOP=../../../..

NMESH=$TOP/exe/nmesh
FDIFF=$TOP/utilities/floatdiff

rm -rf Chkpt_EvryIt

$NMESH advDG-FD.par > advDG-FD.out
$NMESH Chkpt_EvryIt.par > Chkpt_EvryIt.out-1
$NMESH Chkpt_EvryIt.par > Chkpt_EvryIt.out0
$NMESH Chkpt_EvryIt.par > Chkpt_EvryIt.out1
$NMESH Chkpt_EvryIt.par > Chkpt_EvryIt.out2

rm -f Chkpt_AtEnd_correct/Chkpt_AtEnd.par
rm -f Chkpt_EvryIt/Chkpt_EvryIt.par

$FDIFF -at 1e-13 -rt 1e-12 advDG-FD_correct advDG-FD
$FDIFF -at 1e-13 -rt 1e-12 Chkpt_AtEnd_correct Chkpt_EvryIt
# we compare the result of Chkpt_AtEnd.par with Chkpt_EvryIt.par
