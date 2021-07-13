echo "Compiling mcmc_simulation"
mpicc mcmc_simulation.c -o mcmc
echo "Moving executable to scract"
mv mcmc /scratch/castro/local/bin/
