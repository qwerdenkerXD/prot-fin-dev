# bash slurm.sh ../../../materials/protein.fa ../../../materials/mapmanreferencebins.results.txt

function sbatch_script() {
    echo "#!/bin/bash -l

# name
#SBATCH --job-name=ampl_WINSIZE_$3_NPEAKS_$4_AMPL_$5

# cpu
#SBATCH --ntasks=1

#SBATCH --mem-per-cpu=10GB

#SBATCH --output=../results/${exp}/_logs/%x_%j_slurm.out
#SBATCH --error=../results/${exp}/_logs/%x_%j_slurm.err

declare -a overlaps=()
declare -i overlap
for i in 1 2 4; do
    overlap=$3-$3/\${i}
    overlaps+=(\${overlap})
done
overlap=$3-3*$3/4
overlaps+=(\${overlap})
overlap=$3-1
overlaps+=(\${overlap})
for i in \"\${overlaps[@]}\"; do
    name=WINSIZE_$3_NPEAKS_$4_AMPL_$5_OVERLAP_\${i}
    BITS=$5 WINDOW_SIZE=$3 OVERLAP=\${i} N_PEAKS=$4 python3 protfin.py create-db -p ../results/${exp}/\$name.pickle $1
    BITS=$5 WINDOW_SIZE=$3 OVERLAP=\${i} N_PEAKS=$4 python3 protfin.py find-matches -d ../results/${exp}/\$name.pickle ../results/${exp}/_test_selection.fa > ../results/${exp}/_\${name}_test_selection.matches
    awk -v protfin_out=../results/${exp}/_\${name}_test_selection.matches -f extend_protfin_out.awk $2 > ../results/${exp}/_\${name}_test_selection.matches.extended
    rm ../results/${exp}/_\${name}_test_selection.matches
    python3 evaluation.py eval ../results/${exp}/_\${name}_test_selection.matches.extended > ../results/${exp}/_\${name}_test_selection.summary.csv
    python3 evaluation.py plot-extended-out ../results/${exp}/_\${name}_test_selection.matches.extended ../results/${exp}/\$name.png
done
"
}

exp=_v0.3-exp-hashed_amplitudes
mkdir ../results/${exp}
mkdir ../results/${exp}/_logs
python3 evaluation.py select-samples $2 $1 -s 7 > ../results/${exp}/_test_selection.fa
for (( window_size=10; window_size <= 50; window_size+=10 )); do
    for peaks in 5 3 0; do
        for ampl in 1 2 3; do
            sbatch --chdir . <(sbatch_script $1 $2 $window_size $peaks $ampl)
        done
    done
done
