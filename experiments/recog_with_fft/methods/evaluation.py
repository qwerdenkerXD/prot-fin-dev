"""
This module is to do some statistic analysis on the protfin results
"""

from actions import *
import argparse
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL)  # fixes weird python error, look: https://newbebweb.blogspot.com/2012/02/python-head-ioerror-errno-32-broken.html


def main():
    """
    This is where the command line interface is defined to interpret
    the passed arguments to execute the respective functions.
    """

    cli = get_cli()

    args = cli.parse_args()
    args.func(args)


def get_cli():
    parser = argparse.ArgumentParser(
        prog='evaluation',
        description=__doc__,
    )
    parser.set_defaults(func=lambda _: parser.print_help())
    sub_commands = parser.add_subparsers()

    # evaluation.py eval <protfin-out-file>
    eval_parser = sub_commands.add_parser("eval", help="Evaluate protfin find-matches output")
    eval_parser.add_argument("protfin-out-file")
    eval_parser.add_argument("mapman-file")
    eval_parser.set_defaults(func=lambda args: evaluate_protfin(getattr(args, "protfin-out-file"), getattr(args, "mapman-file")))

    # evaluation.py select-samples [-s] <mapman-file> <protein-file>
    eval_parser = sub_commands.add_parser("select-samples", help="Select samples from reference")
    eval_parser.add_argument("mapman-file")
    eval_parser.add_argument("protein-file")
    eval_parser.add_argument("-s", "--samples-per-family", default=1, type=int)
    eval_parser.set_defaults(func=lambda args:
                             select_samples(
                                getattr(args, "mapman-file"),
                                getattr(args, "protein-file"),
                                args.samples_per_family
                             ))

    # evaluation.py plot-extended-out <ext-out-file> <plot-out-file>
    eval_parser = sub_commands.add_parser("plot-extended-out", help="Create a plot of the family extended protfin output")
    eval_parser.add_argument("ext-out-file")
    eval_parser.add_argument("plot-out-file")
    eval_parser.set_defaults(func=lambda args: plot_extended_out(getattr(args, "ext-out-file"), getattr(args, "plot-out-file")))

    # evaluation.py print-hash-counts <database-file>
    eval_parser = sub_commands.add_parser("print-hash-counts", help="Print the calculated hash counts per sequence as csv")
    eval_parser.add_argument("database")
    eval_parser.set_defaults(func=lambda args: print_hash_counts(args.database))

    # evaluation.py print-prots-per-hash <database-file>
    eval_parser = sub_commands.add_parser("print-prots-per-hash", help="Print the counts of proteins for each hash in database as csv")
    eval_parser.add_argument("database")
    eval_parser.set_defaults(func=lambda args: print_prots_per_hash(args.database))

    # evaluation.py plot-frequencies <protein-file> <out-file>
    eval_parser = sub_commands.add_parser("plot-frequencies", help="Plot the selected frequencies")
    eval_parser.add_argument("protein-file")
    eval_parser.add_argument("out-file")
    eval_parser.add_argument("-c", "--cpu", default=1, type=int)
    eval_parser.set_defaults(func=lambda args: plot_frequencies(getattr(args, "protein-file"), getattr(args, "out-file"), args.cpu))

    # evaluation.py plot-prots-per-windist <database-file> <out-file>
    eval_parser = sub_commands.add_parser("plot-prots-per-windist", help="Plot the protein counts per window distance included in the hashes")
    eval_parser.add_argument("database-file")
    eval_parser.add_argument("out-file")
    eval_parser.set_defaults(func=lambda args: plot_prots_per_windist(getattr(args, "database-file"), getattr(args, "out-file")))

    # evaluation.py plot-hashes-per-sequence-length <database-file> <out-file>
    eval_parser = sub_commands.add_parser("plot-hashes-per-sequence-length", help="Plot the hash counts by sequence length")
    eval_parser.add_argument("database-file")
    eval_parser.add_argument("out-file")
    eval_parser.set_defaults(func=lambda args: plot_hashes_per_sequence_length(getattr(args, "database-file"), getattr(args, "out-file")))

    # evaluation.py plot-family-covering <database-file> <out-file>
    eval_parser = sub_commands.add_parser("plot-family-covering", help="Plot the procentual covering of families by their member hashes")
    eval_parser.add_argument("database-file")
    eval_parser.add_argument("mapman-file")
    eval_parser.add_argument("out-file")
    eval_parser.set_defaults(func=lambda args: plot_family_covering(getattr(args, "database-file"), getattr(args, "mapman-file"), getattr(args, "out-file")))

    # evaluation.py plot-hash-frequencies <protein-file> <out-file>
    eval_parser = sub_commands.add_parser("plot-hash-frequencies", help="Plot the average frequencies of hashes in a sequence")
    eval_parser.add_argument("protein-file")
    eval_parser.add_argument("out-file")
    eval_parser.add_argument("-c", "--cpu", default=1, type=int)
    eval_parser.set_defaults(func=lambda args: plot_hash_frequencies(getattr(args, "protein-file"), getattr(args, "out-file"), args.cpu))

    return parser


if __name__ == '__main__':
    main()
