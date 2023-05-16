TOP=../../../..

rm -rf *_correct *_previous

NMESH=$TOP/exe/nmesh
$NMESH sw7FV-DG.par		> sw7FV-DG.out

mv -v sw7FV-DG			sw7FV-DG_correct
