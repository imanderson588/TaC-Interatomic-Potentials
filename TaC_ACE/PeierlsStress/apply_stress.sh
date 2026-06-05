rm energy.yaml
rm -f peierlsstress.out

for i in    80 100 120 140 160 180 200 225 250 275 300 325 350 375 400 425 450 475 500 525 550 ;
do

   mpirun -np 64 lmp -in applystress.in -var stress $i >> peierlsstress.out

done
