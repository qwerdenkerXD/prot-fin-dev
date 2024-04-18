from tools import *
from .algorithm import hashes_from_seq
from multiprocessing import Pool
import pickle


def create_db(
        prot_file: str,
        db_out: str,
        cpu_count=1
        ):
    """
    Creates a database of combinatorial hashes for all protein sequences in the
    given FASTA file and also a lookup database of protein identifiers pointing
    to their description and number of hashes

    ...

    Parameters
    ----------
    prot_file : str
        The path to the FASTA formatted file containing all train sequences of
        known protein sequences
    db_out : str
        Name of the file to write the databases to
    """

    fasta = Fasta(prot_file)

    with Pool(cpu_count) as p:
        subproc_slices = divide_evenly(len(fasta), cpu_count)
        main_slice = next(subproc_slices)
        subprocesses = p.map_async(_process, ((fasta, slc) for slc in subproc_slices))
        database, protein_lookup = _process((fasta, main_slice))

        for sub_db, sub_lookup in subprocesses.get():
            protein_lookup = {**protein_lookup, **sub_lookup}

            for hash_, occs in sub_db.items():
                if hash_ not in database:
                    database[hash_] = occs
                else:
                    database[hash_].extend(occs)

    # write the databases into files
    with open(db_out, 'wb') as db:
        pickle.dump((database, protein_lookup), db, pickle.HIGHEST_PROTOCOL)


def _process(args) ->  Tuple[Database, ProteinLookup]:
    fasta, slc = args

    database: Database = {}
    protein_lookup: ProteinLookup = {}

    # train the database with hashes pointing to their matching proteins
    for prot_id, _, seq in fasta[slc]:

        # calculate the combinatorial hashes for the sequence
        hashes: Hashes = hashes_from_seq(seq, prot_id)

        # save the protein related data
        protein_lookup[prot_id] = (0, len(hashes))

        # sort the hashes into the database by using them as key pointing to
        # their matching proteins
        for hash_, occ in hashes.items():
            if hash_ not in database:
                database[hash_] = []
            database[hash_].append(occ)

    return database, protein_lookup
