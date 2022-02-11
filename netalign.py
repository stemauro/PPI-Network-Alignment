import os
import ppinet as ppin

def ebatch(Ls, taxid):
    getID = lambda x: int(x.split(".")[0])
    return filter(lambda x: getID(x[1]) == taxid, Ls)

def main():

    # Set source directory
    INDIR  = "processed"

    # Read interactions
    fname = "test-virus-host-interactions.txt"
    fpath = os.path.join(os.getcwd(), INDIR, fname)
    Ps = []
    with open(fpath, "r", encoding="utf8") as fin:
        for line in fin:
            tmp = line.rstrip().split("\t")
            Ps.append((tmp[0],tmp[1]))

    # Read similarities
    fname = "test-sim-lcc.txt"
    fpath = os.path.join(os.getcwd(), INDIR, fname)
    Ds = {}
    with open(fpath, "r", encoding="utf8") as fin:
        for line in fin:
            tmp = line.rstrip().split("\t")
            Ds[(tmp[0], tmp[1])] = float(tmp[2])


    # Create networks
    TAXID   = [("10390","10623"),("10390", "10623"),("10390", "10623")]
    COEFF   = [0, 0.5, 1]
    SOLVERS = ["cplex", "glpk"]

    with open("perf.txt", "w") as fout:
        for t in TAXID:
            e1 = ebatch(Ps, int(t[0]))
            e2 = ebatch(Ps, int(t[1]))
            G1 = ppin.PPINet(e1, name=t[0])
            G2 = ppin.PPINet(e2, name=t[1])
            for l in COEFF:
                # Allignment
                alligner = ppin.ILPAlligner(alpha=l)
                alligner.fit(G1, G2, Ds)
                for s in SOLVERS:
                    alligner.solve(solver="glpk", verbose=False)

                    out = "\t".join([alligner.getUID(),
                                min(G1,G2).name, 
                                max(G1,G2).name, 
                                str(alligner.alpha),
                                str(round(alligner.node_similarity(), 4)), 
                                str(round(alligner.edge_correctness(),4)),
                                s,
                                str(alligner.get_wtime())
                               ]) + "\n"
                    fout.write(out)

    pass



if __name__ == "__main__":
    main()
