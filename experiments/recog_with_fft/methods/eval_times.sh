cd "../results/_v0.4-exp-target_zone/_logs"
grep -Po '41669/41669 \[\K[\d:]+' protfin_* |  # extract time from progress bars
sed -E 's/_[0-9]+_slurm.err:|_NPEAKS_|_OVERLAP_|_DIFF_/,/g; s/protfin_WINSIZE_//g' |  # convert to csv: Winsize,n_Peaks,Time
uniq |  # delete duplicates, have the same time, tqdm is not that exact in counting
awk '
BEGIN {
    FS=","
    OFS=","
    prev=""
    print "Winsize,n_Peaks,Overlap,Target_Zone,Time,DB_Size"
}
{
    if ($1 " " $2!=prev) {
        prev=$1 " " $2
        ol = 1
    }
    winsize = $1
    npeaks = $2
    overlap = $3
    diff = $4
    $4 = 2 ** diff
    "du ../protfin_WINSIZE_"winsize"_NPEAKS_"npeaks"_OVERLAP_"overlap"_DIFF_"diff".pickle -hBM | grep -Po \"^\\d+M\"" | getline $6
    print
    ol+=1
}
'
