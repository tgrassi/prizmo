

def preprocess(fname, dictionary):

    in_pragma = False
    fh = open(fname)
    out = ""
    for row in fh:
        srow = row.strip()
        if srow.startswith("!! PREPROCESS_") and srow != "!! PREPROCESS_END":
            nspace = row.split("!!")[0].count(" ")
            if srow.replace("!! PREPROCESS_", "") not in dictionary:
                out += row
                continue
            in_pragma = True
            pragma = dictionary[srow.replace("!! PREPROCESS_", "")]
            indent = " " * nspace
            this_pragma = row + "\n" + indent + pragma.replace("\n", "\n" + indent).rstrip() + "\n\n"
            # this removes the indent in the fortran pragmas (fpp)
            out += "\n".join([x.lstrip() if is_fpp(x) else x for x in this_pragma.split("\n")])
            continue

        if srow == "!! PREPROCESS_END":
            in_pragma = False

        if in_pragma:
            continue

        out += row
    fh.close()

    fout = open(fname, "w")
    fout.write(out)
    fout.close()


# check if the line is a fortran pragma
def is_fpp(line):
    fpps = ["#ifdef", "#else", "#endif"]
    return any([line.strip().startswith(x) for x in fpps])
