exp=_v0.4-exp-uniref_sampling

function sbatch_script() {
    echo "#!/bin/bash -l

# name
#SBATCH --job-name=hashplots_uniref_sample$(basename $pickle)

# cpu
#SBATCH --ntasks=1

#SBATCH --mem-per-cpu=15GB

#SBATCH --output=../results/${exp}/_logs/%x_%j_slurm.out
#SBATCH --error=../results/${exp}/_logs/%x_%j_slurm.err

python3 evaluation.py plot-prots-per-windist \"$pickle\" \"${pickle%.*}_prots_per_windist.png\"
WINDOW_SIZE=${winsize} OVERLAP=${overlap} python3 evaluation.py plot-hashes-per-sequence-length \"$pickle\" \"${pickle%.*}_hashes_per_sequence_length.png\"
python3 evaluation.py plot-family-covering \"$pickle\" ../../../materials/mapmanreferencebins.results.txt \"${pickle%.*}_family_covering.png\" > \"${pickle%.*}_family_covering.csv\"
"
}

for pickle in ../results/${exp}/protfin_*.pickle
do
    winsize=${pickle#*_WINSIZE_}
    winsize=${winsize%%_*}
    overlap=${pickle#*_OVERLAP_}
    overlap=${overlap%%.*}
    sbatch --chdir . <(sbatch_script)
done
