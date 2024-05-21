# bash hashplots_slurm.sh ../../../materials/protein.fa

function sbatch_script() {
    echo "#!/bin/bash -l

# name
#SBATCH --job-name=protfin_WINSIZE_${window_size}_NPEAKS_${peaks}_OVERLAP_${overlap}

# cpu
#SBATCH --ntasks=1

#SBATCH --mem-per-cpu=10GB

#SBATCH --output=../results/${exp}/_logs/%x_%j_slurm.out
#SBATCH --error=../results/${exp}/_logs/%x_%j_slurm.err

db=../results/${exp}/sample_WINSIZE_$window_size.pickle

TITLE=\"Distribution of sequences' hash counts\" X_LABEL=\"Hash counts\" \
Rscript raincloud_plot.R normal <(python3 evaluation.py print-hash-counts \$db) ../results/hash_count_dist.png

python3 evaluation.py plot-frequencies -c 6 $1 ../results/frequencies.png

TITLE=\"Distribution of proteins per hash\" X_LABEL=\"Protein counts\" \
Rscript raincloud_plot_log10.R normal <(python3 evaluation.py print-prots-per-hash \$db) ../results/prots_per_hash.png

python3 evaluation.py plot-prots-per-windist \$db ../results/prots_per_windist.png
"
}

exp=_v0.4-exp-uniref_sampling

for (( window_size=10; window_size <= 50; window_size+=10 )); do
    declare -a overlaps=()
    declare -i overlap
    for i in 1 2 4; do
        overlap=${window_size}-${window_size}/${i}
        overlaps+=(${overlap})
    done
    overlap=${window_size}-3*${window_size}/4
    overlaps+=(${overlap})
    overlap=${window_size}-1
    overlaps+=(${overlap})

    for overlap in "${overlaps[@]}"; do
        for peaks in 5 3 0; do
            sbatch --chdir . <(sbatch_script $*)
        done
    done
done
