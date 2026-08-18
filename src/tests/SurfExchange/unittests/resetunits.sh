TOP=../../../..

rm -rf *_correct *_previous

NMESH=$TOP/exe/nmesh

/bin/cp SurfExchange_unit1.par SurfExchange_unit1_p1.par
/bin/cp SurfExchange_unit1.par SurfExchange_unit1_p2.par
/bin/cp SurfExchange_unit1.par SurfExchange_unit1_p3.par

rm -rf SurfExchange_unit1_p?_corr

mpirun -np 1 $NMESH SurfExchange_unit1_p1.par > SurfExchange_unit1_p1.out
mv -v SurfExchange_unit1_p1 SurfExchange_unit1_p1_corr

mpirun -np 2 $NMESH SurfExchange_unit1_p2.par > SurfExchange_unit1_p2.out
mv -v SurfExchange_unit1_p2 SurfExchange_unit1_p2_corr

mpirun -np 3 $NMESH SurfExchange_unit1_p3.par > SurfExchange_unit1_p3.out
mv -v SurfExchange_unit1_p3 SurfExchange_unit1_p3_corr
