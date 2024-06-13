# bash protfin_slurm.sh ../../../materials/protein.fa ../../../materials/mapmanreferencebins.results.txt

function sbatch_script() {
    name=../results/${exp}/protfin_WINSIZE_${window_size}_NPEAKS_${peaks}_OVERLAP_${overlap}
    filtname=${name}_FILT_${filt}
    echo "#!/bin/bash -l

# name
#SBATCH --job-name=protfin_WINSIZE_${window_size}_NPEAKS_${peaks}_OVERLAP_${overlap}_FILT_${filt}

# cpu
#SBATCH --ntasks=1

#SBATCH --mem-per-cpu=15GB

#SBATCH --output=../results/${exp}/_logs/%x_%j_slurm.out
#SBATCH --error=../results/${exp}/_logs/%x_%j_slurm.err

WINDOW_SIZE=${window_size} OVERLAP=${overlap} N_PEAKS=${peaks} python3 protfin.py create-db -p ${name}.pickle $1
WINDOW_SIZE=${window_size} OVERLAP=${overlap} N_PEAKS=${peaks} python3 protfin.py find-matches -f ${filt} -d ${name}.pickle ../results/${exp}/_test_selection.fa > ${filtname}.matches
awk -v protfin_out=${filtname}.matches -f extend_protfin_out.awk $2 > ${filtname}.matches.extended
rm ${filtname}.matches
python3 evaluation.py eval ${filtname}.matches.extended > ${filtname}.summary.csv
"
}

exp=_v0.4-exp-filter_hashes
mkdir ../results/${exp}
mkdir ../results/${exp}/_logs
python3 evaluation.py select-samples $2 $1 -s 7 > ../results/${exp}/_test_selection.fa

for (( window_size=50; window_size >= 10; window_size-=10 )); do
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
            for filt in .1 .25 .5 .75; do
                sbatch --chdir . <(sbatch_script $*)
            done
        done
    done
done
