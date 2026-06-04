rm energy.yaml

for i in    80 100 120 140 160 180 200 225 250 275 300 ;
do

   mpirun -np 64 lmp -in applystress.in -var stress $i

done
