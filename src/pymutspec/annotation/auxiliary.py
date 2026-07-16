transcriptor = str.maketrans("ACGT", "TGCA")


def rev_comp(mut: str):
    """
    Return the reverse complement of a 192-component SBS mutation string.

    The input format is ``X[N>M]Y`` where ``X`` and ``Y`` are the flanking
    nucleotides and ``N>M`` is the substitution.  The function swaps the
    flanking nucleotides and applies complement translation to all characters.

    Arguments
    ---------
    mut: str
        SBS mutation string of the form ``X[N>M]Y``.

    Return
    ------
    rev_comp_mut: str
        Reverse-complemented mutation string.
    """
    new_mut = mut[-1] + mut[1:-1] + mut[0]
    new_mut = new_mut.translate(transcriptor)
    return new_mut


def lbl_id2lbl(lbl_id: int) -> str:
    """
    Convert a numeric label identifier to its string label.

    Arguments
    ---------
    lbl_id: int
        Integer label code: 0 → ``'all'``, 1 → ``'syn'``, 2 → ``'ff'``.

    Return
    ------
    lbl: str
        Human-readable label string.

    Raises
    ------
    NotImplementedError
        If ``lbl_id`` is not 0, 1, or 2.
    """
    if lbl_id == 0:
        lbl = "all"
    elif lbl_id == 1:
        lbl = "syn"
    elif lbl_id == 2:
        lbl = "ff"
    else:
        raise NotImplementedError()
    return lbl


def lbl2lbl_id(lbl: str) -> int:
    """
    Convert a string label to its numeric identifier.

    Arguments
    ---------
    lbl: str
        Label string; one of ``'all'``, ``'syn'``, ``'syn_c'``,
        ``'ff'``, or ``'syn4f'``.

    Return
    ------
    lbl_id: int
        Integer label code: ``'all'`` → 0, ``'syn'``/``'syn_c'`` → 1,
        ``'ff'``/``'syn4f'`` → 2.

    Raises
    ------
    NotImplementedError
        If ``lbl`` is not a recognised label string.
    """
    if lbl == "all":
        lbl_id = 0
    elif lbl == "syn" or lbl == "syn_c":
        lbl_id = 1
    elif lbl == "ff" or lbl == "syn4f":
        lbl_id = 2
    else:
        raise NotImplementedError()
    return lbl_id
