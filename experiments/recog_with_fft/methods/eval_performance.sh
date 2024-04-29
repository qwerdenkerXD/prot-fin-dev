cd ../results/*/_logs
grep -Po '41669/41669 \[\K[\d:]+' peaksel_* |  # extract time from progress bars
sed -E 's/_[0-9]+_slurm.err:|_NPEAKS_/,/g; s/peaksel_WINSIZE_//g' |  # convert to csv: Winsize,n_Peaks,Time
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
    time = $3
    overlaps[1] = 0
    overlaps[2] = winsize - int(winsize/2)
    overlaps[3] = winsize - int(winsize/4)
    overlaps[4] = winsize - int(3*winsize/4)
    overlaps[5] = winsize - 1
    overlap = overlaps[ol]
    $3=overlap
    $4=time
    "du ../WINSIZE_"winsize"_NPEAKS_"npeaks"_OVERLAP_"overlap".pickle -hBM | grep -Po \"^\\d+M\"" | getline $5
    print
    ol+=1
}
'