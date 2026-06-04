rm gsf_111.yaml
rm gsf_112.yaml 
rm gsf_110.yaml 


for i in 0.000 0.05 0.1 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.0 #1.05 1.10 1.15 1.20 1.25 1.30 1.35 1.40 1.45 1.50 1.55 1.60 1.65 1.70 1.75 1.80 1.85 1.90 1.95 2.0 2.2 2.3 2.4 2.5 2.6
do
   echo "Running $i"
   
   mpirun -np 64 lmp -in gsf_112.in -var x_shift $i
   mpirun -np 64 lmp -in gsf_110.in -var x_shift $i
   mpirun -np 64 lmp -in gsf_111.in -var x_shift $i

done
