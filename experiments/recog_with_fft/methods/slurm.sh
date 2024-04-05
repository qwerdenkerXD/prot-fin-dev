function sbatch_script() {
    echo "#!/bin/bash -l

# name
#SBATCH --job-name=WINSIZE_$3_OVERLAP_$4_NPEAKS_$5

# cpu
#SBATCH --ntasks=1

#SBATCH --mem-per-cpu=1GB

#SBATCH --output=../results/_stft_param_exp/logs/%x_%j_slurm.out
#SBATCH --error=../results/_stft_param_exp/logs/%x_%j_slurm.err

name=WINSIZE_$3_OVERLAP_$4_NPEAKS_$5
WINDOW_SIZE=$3 OVERLAP=$4 N_PEAKS=$5 python3 protfin.py create-db -p ../results/_stft_param_exp/\$name.pickle $1
WINDOW_SIZE=$3 OVERLAP=$4 N_PEAKS=$5 python3 protfin.py find-matches -d ../results/_stft_param_exp/\$name.pickle ../results/_stft_param_exp/_test_selection.fa > ../results/_stft_param_exp/_\${name}_test_selection.matches
python3 evaluation.py eval ../results/_stft_param_exp/_\${name}_test_selection.matches > ../results/_stft_param_exp/_\${name}_test_selection.summary.csv
"
}

mkdir ../results/_stft_param_exp
mkdir ../results/_stft_param_exp/logs
python3 evaluation.py select-samples $2 $1 -s 7 > ../results/_stft_param_exp/_test_selection.fa
for (( window_size=10; window_size <= 50; window_size+=5 )); do
    for hop_size in 2 4 5 10; do
        declare -i overlap=$window_size-$window_size/$hop_size
        for peaks in 5 3 0; do
            sbatch --chdir . <(sbatch_script $1 $2 $window_size $overlap $peaks)
        done
    done
done
