exp=_v0.4-exp-uniref_sampling

function sbatch_script() {
    echo "#!/bin/bash -l

# name
#SBATCH --job-name=plot_uniref_sample_ca_150Mio

# cpu
#SBATCH --ntasks=1

#SBATCH --mem-per-cpu=15GB

#SBATCH --output=../results/${exp}/_logs/%x_%j_slurm.out
#SBATCH --error=../results/${exp}/_logs/%x_%j_slurm.err

for pickle in ../results/${exp}/sample_*.pickle
do
    while :
    do
        python3 uniref_sample_exp.py plot \"\$pickle\" \"\${pickle%.*}.png\" > \"\${pickle%.*}.csv\"&& break  # can be that pickle file changes during plotting -> error
    done
done
"
}

sbatch --chdir . <(sbatch_script)
