import os
import ppinet as ppin
from itertools import combinations

def ebatch(Ls, taxid):
    getID = lambda x: int(x.split(".")[0])
    return filter(lambda x: getID(x[1]) == taxid, Ls)

def main():

    # Set source directory
    INDIR  = "processed"

    # Read interactions
    fname = "exam-virus-host-interactions.txt"
    fpath = os.path.join(os.getcwd(), INDIR, fname)
    Ps = []
    with open(fpath, "r", encoding="utf8") as fin:
        for line in fin:
            tmp = line.rstrip().split("\t")
            Ps.append((tmp[0],tmp[1]))

    # Read similarities
    fname = "exam-sim-lcc.txt"
    fpath = os.path.join(os.getcwd(), INDIR, fname)
    Ds = {}
    with open(fpath, "r", encoding="utf8") as fin:
        for line in fin:
            tmp = line.rstrip().split("\t")
            Ds[(tmp[0], tmp[1])] = float(tmp[2])


    # Create networks
    TAXID   = ["11552", "334203", "11070", "11082", "10255", "10407", "138950", "11520"]
    SOLVERS = ["cplex", "glpk"]

    #t = TAXID[0]
    with open("perf.txt", "w") as fout:
        for t in combinations(TAXID,2):
            coeff = 0
            e1 = ebatch(Ps, int(t[0]))
            e2 = ebatch(Ps, int(t[1]))
            G1 = ppin.PPINet(e1, name=t[0])
            G2 = ppin.PPINet(e2, name=t[1])
            alligner = ppin.ILPAlligner()
            alligner.fit(G1, G2, Ds, coeff)
            while coeff  <= 1:
                alligner.update_obj(coeff)
                for s in SOLVERS:
                    alligner.solve(solver=s, verbose=False)
                    out = "\t".join([alligner.getUID(),
                                min(G1,G2).name, 
                                max(G1,G2).name, 
                                str(coeff),
                                str(round(alligner.node_similarity(), 4)), 
                                str(round(alligner.edge_correctness(),4)),
                                s,
                                str(alligner.get_wtime())
                               ]) + "\n"
                    fout.write(out)
                coeff += 0.25

    pass



if __name__ == "__main__":
    main()
