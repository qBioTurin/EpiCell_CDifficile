
# usage for analysis:
# * navigate on Main.sh's dir (./Global_SA)
# * [bash]: chmod +x Main.sh
# * [bash]: ./Main.sh

wd="/home/raucello/EpiCell_CDifficile"
cores=1
node=4
samples=$(($cores*$node))

model_mat="Recon3D.mat"

cd $wd/Notebooks/Global_SA/
mpiexec -n $cores python3 saltelli_sample.py ${model_mat}  --num_samples $samples
python3 sobol_analyze.py

python3 $wd/code/global_sens_analysis/pickle_reader.py