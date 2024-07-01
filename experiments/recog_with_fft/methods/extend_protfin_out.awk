# awk -v protfin_out=<protfin-out.csv> -f extend_protfin_out.awk <mapman-file>

BEGIN {
    FS=OFS="\t"
}
{
    if ($3 != "IDENTIFIER" && $3 != "''") {
        prot_id = substr(toupper($3), 2, length($3)-2)
        bin = substr(toupper($1), 2, length($1)-2)
        if (bins[prot_id]) {
            bins[prot_id] = bins[prot_id]"|"bin
        } else {
            bins[prot_id] = bin
        }
    }
}
END {
    if (!protfin_out) {
        print "Missing variable, please set with -v option:"
        print "awk -v protfin_out=<protfin-out.csv> ..."
        exit
    }
    while ((getline line < protfin_out) > 0) {
        len=split(line, fields, ",")
        if (fields[1] == "Rank") {
            print line",Input_Family,Match_Family"
        } else {
            match_id = fields[2]
            input_id = fields[5]
            print line","bins[input_id]","bins[match_id]
        }
    }
    close(protfin_out)
}
