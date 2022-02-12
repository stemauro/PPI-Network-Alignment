import os
from itertools import groupby

def ishh(l, taxid):
    idfir, idsec = int(l[0].split(".")[0]), int(l[1].split(".")[0])
    return idfir == idsec == taxid

def isbelow(l, t):
    score = int(l[-1])
    return score < t
    
def isin(l, Ls):
    idref = int(l[1].split(".")[0])
    return idref in Ls

def protid(l):
    return l.split(" ")[0][1:]

def main() :

    # Environment variables
    INDIR               = "files"
    OUTDIR              = "processed" 
    DESTDIR             = os.path.join(os.getcwd(), OUTDIR)
    THRESHOLD           = 510
    HOMO_SAPIENS_TAXID  = 9606
    SOL                 = ">"
    #PASS_IDS            = [11269,186538,194443,194441,11137,194440,
    #                       694009,337041,11320,11103,162145,11250,
    #                       32604,32603,11161,63330,11234,11676,11709,
    #                       37296,10359,10376,10335,10310,10298]
    #PASS_IDS            = [10390, 10623]
    PASS_IDS            = [11552, 334203, 11070, 11082,
                           10255, 10407, 138950, 11520]

    # I/O filepaths
    infile  = "9606.protein.links.v10.5.txt"
    outfile = "exam-virus-host-interactions.txt"
    inpath  = os.path.join(os.getcwd(), INDIR, infile)
    outpath = os.path.join(os.getcwd(), OUTDIR, outfile)

    # Create target directory if non-existent
    if not os.path.isdir(DESTDIR):
        os.mkdir(DESTDIR)
     
    # Parse virus-host protein-protein interactions
    fin  = open(inpath, "r", encoding="utf8")
    fout = open(outpath, "w")
    Ps = []

    print("Parsing interactions from", inpath)
    fin.readline()  # skip first row
    for line in fin:
        tmp = line.split(" ")
        if ishh(tmp, HOMO_SAPIENS_TAXID) or isbelow(tmp, THRESHOLD) or not isin(tmp, PASS_IDS):
            continue
        Ps.append(tmp[0])
        Ps.append(tmp[1])
        fout.write(tmp[0] + "\t" + tmp[1] + "\n")
    fout.close()
    fin.close()
  
    # Unique proteins extracted
    pass_prot_ids = [k for k, _ in groupby(sorted(Ps))]
    print("Unique proteins extracted:", len(pass_prot_ids))
    print("Output file:", outpath)
    
    #raise KeyboardInterrupt()
    
    # I/O filepaths
    infile  = "protein.sequences.v10.5.fa"
    outfile = "exam-protein-sequences.fa"
    inpath  = os.path.join(os.getcwd(), INDIR, infile)
    outpath = os.path.join(os.getcwd(), OUTDIR, outfile)
    
    # Parse aminoacid sequences
    fin  = open(inpath, "r", encoding="utf8")
    fout = open(outpath, "w")
    wflag = False
    count = 0

    print("Parsing sequences from", inpath)
    for line in fin:
        if line[0] == SOL:
            wflag = protid(line) in pass_prot_ids
            count += wflag
        if wflag:
            fout.write(SOL + protid(line) + "\n") if line[0] == SOL else fout.write(line)
    fout.close()
    fin.close()
    print("Number of parsed sequences:", count)
    print("Output file:", outpath)

    pass



if __name__ == "__main__":

    main()
