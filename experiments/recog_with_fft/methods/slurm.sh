# bash slurm.sh ../../../materials/protein.fa

function sbatch_script() {
    echo "#!/bin/bash -l

# name
#SBATCH --job-name=sample_WINSIZE_$2

# cpu
#SBATCH --ntasks=1

#SBATCH --mem-per-cpu=10GB

#SBATCH --output=../results/${exp}/_logs/%x_%j_slurm.out
#SBATCH --error=../results/${exp}/_logs/%x_%j_slurm.err

WINDOW_SIZE=$2 python3 uniref_sample_exp.py run $1 ../results/${exp}/sample_WINSIZE_$2.pickle
"
}

exp=_v0.4-exp-uniref_sampling
mkdir ../results/${exp}
mkdir ../results/${exp}/_logs
for (( window_size=10; window_size <= 100; window_size+=10 )); do
    sbatch --chdir . <(sbatch_script $1 $window_size)
done
