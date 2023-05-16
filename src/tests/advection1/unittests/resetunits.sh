TOP=../../../..

rm -rf *_correct *_previous

NMESH=$TOP/exe/nmesh
$NMESH advDG-FD.par		> advDG-FD.out

mv -v advDG-FD			advDG-FD_correct
