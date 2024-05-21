cd ../results/_v0.4*/_logs
grep -Po '41669/41669 \[\K[\d:]+' protfin_* |  # extract time from progress bars
sed -E 's/_[0-9]+_slurm.err:|_NPEAKS_|_OVERLAP_/,/g; s/protfin_WINSIZE_//g' |  # convert to csv: Winsize,n_Peaks,Time
uniq |  # delete duplicates, have the same time, tqdm is not that exact in counting
awk '
BEGIN {
    FS=","
    OFS=","
    prev=""
    print "Winsize,n_Peaks,Overlap,Time,DB_Size"
}
{
    if ($1 " " $2!=prev) {
        prev=$1 " " $2
        ol = 1
    }
    winsize = $1
    npeaks = $2
    overlap = $3
    "du ../protfin_WINSIZE_"winsize"_NPEAKS_"npeaks"_OVERLAP_"overlap".pickle -hBM | grep -Po \"^\\d+M\"" | getline $5
    print
    ol+=1
}
'
