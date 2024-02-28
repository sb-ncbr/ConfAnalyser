def load_file(file_name: str) -> list[str]:
    try:
        with open(file_name) as file:
            return file.readlines()
    except OSError:
        print(f"[ERROR] File `{file_name}` could not be opened!\nCheck the path and try again.")
        exit()


def load_names(file_name: str) -> dict[str, list[set[str]]]:
    lines = load_file(file_name)
    out: dict[str, list[set[str]]] = dict()
    for line in lines:
        split_line = line.split()
        ligand = split_line[0]
        if ligand[2] == "_":
            ligand = ligand[0:2]
        if ligand in out:
            for i, name in enumerate(split_line[1:]):
                out[ligand][i].add(name)
        else:
            out[ligand] = []
            for name in split_line[1:]:
                out[ligand].append({name})
    return out

def print_dict(dct: dict[str, set[str]]) -> None:
    """
    Prints out a debug output of a dictionary.
    Used to print list of names of atoms loaded from atom_names.txt
    """
    for name in dct:
        print(f"Ligand: {name}")
        for lst in dct[name]:
            out = ", ".join(x for x in lst)
            print(f"   -> {out}")