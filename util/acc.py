def read_accession_file(f):
    """
    Read an open accession file, returning the list of accession numbers it
    contains.

    This automatically skips blank lines and comments.
    """
    for l in f:
        l = l.strip()
        if not l or l.startswith('#'):
            continue
        yield l

