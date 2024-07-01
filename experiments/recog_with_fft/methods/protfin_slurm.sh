# bash protfin_slurm.sh ../../../materials/protein.fa ../../../materials/mapmanreferencebins.results.txt

function sbatch_script() {
    name=protfin_WINSIZE_${window_size}_NPEAKS_${peaks}_OVERLAP_${overlap}_FIRSTAPPEARANCE_${first_appearance}
    result_name=../results/${exp}/${name}
    echo "#!/bin/bash -l

# name
#SBATCH --job-name=${name}

# cpu
#SBATCH --ntasks=1

#SBATCH --mem-per-cpu=15GB

#SBATCH --output=../results/${exp}/_logs/%x_%j_slurm.out
#SBATCH --error=../results/${exp}/_logs/%x_%j_slurm.err

echo createdb >&2
FIRST_APPEARANCE=${first_appearance} WINDOW_SIZE=${window_size} OVERLAP=${overlap} N_PEAKS=${peaks} python3 protfin.py create-db -p ${result_name}.pickle $1

for filt in .1 .2 .3 .4 .5 .6 .7 .8. .9 1; do
    echo findmatches FILT \$filt >&2
    FIRST_APPEARANCE=${first_appearance} WINDOW_SIZE=${window_size} OVERLAP=${overlap} N_PEAKS=${peaks} python3 protfin.py find-matches -f \${filt} -d ${result_name}.pickle ../results/${exp}/_test_selection.fa > ${result_name}_FILT_\${filt}.matches
    echo eval FILT \$filt >&2
    python3 evaluation.py eval ${result_name}_FILT_\${filt}.matches $2 > ${result_name}_FILT_\${filt}.summary.csv
    echo matchfamily FILT \$filt >&2
    python3 protfin.py match-family <(awk \"BEGIN{FS=\\\"\t\\\"; print \\\"Family_ID,Protein_ID\\\"} \\\$5==\\\"T\\\"{x=\\\$1\\\",\\\"toupper(\\\$3);gsub(/'/, \\\"\\\", x);print x}\" $2) -d ${result_name}.pickle -f \${filt} > ${result_name}_FILT_\${filt}.match_fam.csv
done

echo plothashfrequencies >&2
FIRST_APPEARANCE=${first_appearance} python3 evaluation.py plot-hash-frequencies $1 ${result_name}.png > ${result_name}.hash_frequencies.csv
"
}

exp=_v0.4-exp-filter_hashes
mkdir ../results/${exp}
mkdir ../results/${exp}/_logs
python3 evaluation.py select-samples $2 $1 -s 7 > ../results/${exp}/_test_selection.fa

for (( window_size=50; window_size >= 20; window_size-=10 )); do
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
        for peaks in 5 3; do
            for first_appearance in 0 1; do
                sbatch --chdir . <(sbatch_script $*)
            done
        done
    done
done
