def read_accession_file(f):
    """
    Read an open accession file, returning the list of accession numbers it
    contains.

    This automatically skips blank lines and comments.
    """
    for line in f:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        yield line
