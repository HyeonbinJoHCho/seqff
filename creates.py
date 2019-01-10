import os


def create_newtemp(newtemp, file, outdir):
    """
    Returns nothing
    Creates a counts-formatted file as outpath using bincount(=newtemp) dictionary data

    :param bincount: dict
    :param file: filename of sam file
    :param outdir: str
    """
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    outfilename = os.path.join(outdir, os.path.split(file)[1]) + ".newtemp"
    with open(outfilename ,'w') as outfile:
        for key, val in newtemp.items():
            print(key, val, sep="\t", file=outfile)


def create_normalizedbincount(normalizedbincount, file, outdir):
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    outfilename = os.path.join(outdir, os.path.split(file)[1]) + ".normalizedbincount"
    with open(outfilename, 'w') as outfile:
        for data in normalizedbincount:
            print(data, file=outfile)


def create_alluseablebins(alluseablebins, file, outdir):
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    outfilename = os.path.join(outdir, os.path.split(file)[1]) + ".alluseablebins"
    with open(outfilename, 'w') as outfile:
        for data in alluseablebins:
            print(data, file=outfile)


def create_bincounts(bincounts, file, outdir):
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    outfilename = os.path.join(outdir, os.path.split(file)[1]) + ".bincounts"
    with open(outfilename, 'w') as outfile:
        for data in bincounts:
            print(data, file=outfile)


def create_results(lines, file, outdir):
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    outfilename = os.path.join(outdir, os.path.split(file)[1]) + ".results"
    with open(outfilename, 'w') as outfile:
        for line in lines:
            print(line, file=outfile)
