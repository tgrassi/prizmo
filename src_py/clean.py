from glob import glob
import os

cwd = os.getcwd()
if not cwd.endswith("src_py"):
    print("Please run this script from the src_py directory.")
    exit(1)

makefile_pragmas = "../Makefile_pragmas"
if os.path.isfile(makefile_pragmas):
    os.remove(makefile_pragmas)
    print("removed " + makefile_pragmas)

for g in glob("../*f90"):
    print("cleaning " + g)
    s = ""
    in_pragma = False
    fh = open(g)
    for row in fh:
        if row.strip().startswith("!! PREPROCESS_"):
            in_pragma = not in_pragma
            s += row
            continue
        if in_pragma:
            continue
        s += row
    fh.close()

    open(g, "w").write(s)
