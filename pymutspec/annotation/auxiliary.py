transcriptor = str.maketrans("ACGT", "TGCA")


def rev_comp(mut: str):
    new_mut = mut[-1] + mut[1:-1] + mut[0]
    new_mut = new_mut.translate(transcriptor)
    return new_mut


def lbl_id2lbl(lbl_id: int) -> str:
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
    if lbl == "all":
        lbl_id = 0
    elif lbl == "syn" or lbl == "syn_c":
        lbl_id = 1
    elif lbl == "ff" or lbl == "syn4f":
        lbl_id = 2
    else:
        raise NotImplementedError()
    return lbl_id
