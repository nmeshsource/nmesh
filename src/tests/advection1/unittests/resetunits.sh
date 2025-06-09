TOP=../../../..

rm -rf *_correct *_previous  Chkpt_AtEnd

NMESH=$TOP/exe/nmesh
$NMESH advDG-FD.par		> advDG-FD.out
$NMESH Chkpt_AtEnd.par		> Chkpt_AtEnd.out

mv -v advDG-FD			advDG-FD_correct
mv -v Chkpt_AtEnd		Chkpt_AtEnd_correct
# we intend to compare the result of Chkpt_AtEnd.par with Chkpt_EvryIt.par
