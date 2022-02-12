import sys
import random
import string
import networkx as nx
import numpy as np
import pyomo.environ as pyo
from itertools import product
from scipy.sparse import csr_matrix

class PPINet(nx.Graph):

    def __gt__(self, other):
        return len(self.nodes) > len(other.nodes)
    
    pass

def _getsim(G1, G2, Ds):
    Ss = {}
    V1 = list(G1.nodes)
    V2 = list(G2.nodes)
    r1 = range(len(V1))
    r2 = range(len(V2))
    for i, j in product(r1, r2):
        try:
            Ss[(i,j)] = Ds[(V1[i],V2[j])]
        except:
            try:
                Ss[(i,j)] = Ds[(V2[j],V1[i])]
            except:
                if V1[i] == V2[j]:
                    Ss[(i,j)] = 1.0
                else:
                    raise KeyError("({},{})".format(i,j))

    return Ss

class IQPAlligner(object):
    """
    PPI network allignment solved as an integrer quadratic programming (IQP) problem.
    
    The optimization objective for IQPAlligner is:

    

    obtaines as a convex combination of two scores measuring
    biological coherence and topological coherence, respectively.

    Values of functional coherence (FC) and edge correctness (EC)
    are obtained by normalizing biological coherence and
    topological coherence scores by the number of nodes and edges
    of the smallest network, respectively.

    """
    def __init__(self, alpha, time_limit=60, prob_file=False):
        
        # Safe check on parameters
        if alpha < 0 or alpha > 1:
            raise ValueError("The alpha coefficient must be between 0 and 1.")
        self._uid       = ''.join(random.choices(string.ascii_uppercase + string.digits, k=8))
        self._timeout   = time_limit
        self._solvers   = ["cplex", "gurobi", "ipopt"]
        self._prob_file = prob_file 
        self.alpha      = alpha

    def getUID(self):
        return self._uid

    def fit(self, G1, G2, Ds):

        # Ensure smaller networks are mapped into
        # larger networks for consistency
        if G1 > G2:
            G1, G2 = G2, G1
        
        # Compute adjacency matrix of each network
        # and similarity scores between networks
        self.A, self.B = map(nx.adjacency_matrix, [G1, G2])
        self.S =_getsim(G1, G2, Ds)

        # Number of nodes in each network
        self._n1 = len(G1.nodes)
        self._n2 = len(G2.nodes)

        print("Generating problem instance (IID: {})".format(self._uid), file=sys.stderr)

        # Define IQP model
        model = pyo.ConcreteModel()
        model.I = pyo.RangeSet(0, self._n1 - 1)
        model.K = pyo.RangeSet(0, self._n2 - 1)

        # Decision variables
        model.x = pyo.Var(model.I, model.K, domain=pyo.Binary)

        # Objective function
        def bc(model, Ps):
            return sum(Ps[(i,k)] * model.x[i,k] for k in model.K for i in model.I)

        def tc(model, A1, A2):
            return sum(A1[i,j] * A2[k,l] * model.x[i,k] * model.x[j,l]
                    for l in model.K for j in model.I
                    for k in model.K for i in model.I)
        
        def obj_rule(model):
            fc = bc(model, self.S)
            ec = tc(model, self.A, self.B)
            return self.alpha * fc + (1 - self.alpha) * ec

        model.obj = pyo.Objective(rule=obj_rule, sense=-1)
        

        # X is a (partial) mapping
        def cst_mapping(model, i):
            return sum(model.x[i,k] for k in model.K) <= 1
        model.mapping = pyo.Constraint(model.I, rule=cst_mapping)

        # X is injective
        def cst_injective(model, k):
            return sum(model.x[i,k] for i in model.I) <= 1
        model.injective = pyo.Constraint(model.K, rule=cst_injective)
        
        # Write model to file
        if self._prob_file:
            outpath = os.path.abspath(''.join(["iqp",self._uid,".lp"]))
            model.write(outpath)
        
        self.model = model


    def solve(self, solver="ipopt", verbose=True):
        s = solver.lower()
        if s not in self._solvers:
            raise ValueError("Unsupported solver requested. Available solvers are: {}".format(self._solvers))
        
        sol = pyo.SolverFactory(s)
        if self._timeout is not None:
            # If specified, set solver timeout.
            #
            # Default timeout is one minute regardless of the solver being invoked.
            # Unsetting the time limit is highly discouraged, given the
            # computational burden of exact algorithms (like this)
            # even for medium-sized instances of the allignment problem.
            if s == "cplex":
                sol.options["timelimit"] = self._timeout
            elif s == "glpk":
                sol.options["tmlim"] = self._timeout
            elif s == "gurobi":
                sol.options["TimeLimit"] = self._timeout 
        
        print("Solving instance...", file=sys.stderr)

        self.solved_model = sol.solve(self.model, tee=verbose)
        self._sol_json = self.solved_model.json_repn()

        # Check solver status
        status = self._sol_json["Solver"][0]["Status"]
        if status == "ok":
            data = np.array([self.model.x[i,k]() for i in self.model.I for k in self.model.K])
            self.opt = csr_matrix(data.reshape((self._n1, self._n2)))
            self.opt_val = self.model.obj()
        else:
            raise RuntimeError("Solver status:", status)

    def node_similarity(self):
        bc = sum(self.S[(i,k)] * self.model.x[i,k]() for k in self.model.K for i in self.model.I)
 
        return bc / self._n1
    
    def edge_correctness(self):
        tc = sum(self.A[i,j] * self.B[k,l] * self.model.x[i,k]() * self.model.x[j,l]()
                for l in self.model.K for j in self.model.I
                for k in self.model.K for i in self.model.I)
    
        return tc / self.A.sum()

    def get_wtime(self):
        return self.solved_model.solver.wallclock_time


class ILPAlligner(object):
    """
    PPI network allignment solved as an integrer linear programming (ILP) problem.
    
    The optimization objective for ILPAlligner is:
    
    obtaines as a convex combination of two scores measuring
    biological coherence and topological coherence, respectively.
    Values of functional coherence (FC) and edge correctness (EC)
    are obtained by normalizing biological coherence and
    topological coherence scores by the number of nodes and edges
    of the smallest network, respectively.
    """
    def __init__(self, time_limit=60, prob_file=False):
        self._uid       = ''.join(random.choices(string.ascii_uppercase + string.digits, k=8))
        self._timeout   = time_limit
        self._solvers   = ["cplex", "glpk", "gurobi"]
        self._prob_file = prob_file 

    def getUID(self):
        return self._uid

    def _obj_rule(self, m, alpha):
        fc = sum(self.S[(i,k)] * m.x[i,k] for k in m.K for i in m.I)
        ec = sum(m.y[i,k] for k in m.K for i in m.I)
        return alpha * fc + (1 - alpha) * ec


    def fit(self, G1, G2, Ds, alpha):

        if alpha < 0 or alpha > 1:
            raise ValueError("The objective coefficient must be between 0 and 1.")

        # Ensure smaller networks are mapped into
        # larger networks for consistency
        if G1 > G2:
            G1, G2 = G2, G1
        
        # Compute adjacency matrix of each network
        # and similarity scores between networks
        self.A, self.B = map(nx.adjacency_matrix, [G1, G2])
        self.S = _getsim(G1, G2, Ds)

        # Number of nodes in each network
        self._n1 = len(G1.nodes)
        self._n2 = len(G2.nodes)

        print("Generating problem instance (UID: {})".format(self._uid), file=sys.stderr)
        
        # Define IQP model
        model = pyo.ConcreteModel()
        model.I = pyo.RangeSet(0, self._n1 - 1)
        model.K = pyo.RangeSet(0, self._n2 - 1)

        # Decision variables
        model.x = pyo.Var(model.I, model.K, domain=pyo.Binary)
        model.y = pyo.Var(model.I, model.K, domain=pyo.NonNegativeIntegers)

        # X is a (partial) mapping
        def cst_mapping(model, i):
            return sum(model.x[i,k] for k in model.K) <= 1
        model.mapping = pyo.Constraint(model.I, rule=cst_mapping)

        # X is injective
        def cst_injective(model, k):
            return sum(model.x[i,k] for i in model.I) <= 1
        model.injective = pyo.Constraint(model.K, rule=cst_injective)

        # Relationship between x and y
        def cst_bin2int(model, i, k):
            c = sum(self.A[i,j] * self.B[k,l] for l in model.K for j in model.I)
            return model.y[i,k] <= c * model.x[i,k]
        model.bin2int = pyo.Constraint(model.I, model.K, rule=cst_bin2int)

        # Number of preserved edges
        def cst_preserved(model, i, k):
            ub = sum(self.A[i,j] * self.B[k,l] * model.x[j,l] for l in model.K for j in model.I)
            return model.y[i,k] <= ub
        model.preserved = pyo.Constraint(model.I, model.K, rule=cst_preserved)
        
        # Write model to file
        if self._prob_file:
            outpath = os.path.abspath(''.join(["ilp",self._uid,".lp"]))
            model.write(outpath)
       
        # Objective function
        model.obj = pyo.Objective(rule=lambda x: self._obj_rule(x, alpha), sense=pyo.maximize)
        
        self.model = model

        pass
    
    
    def update_obj(self, alpha):
        if hasattr(self, "model"):
            self.model.del_component(self.model.obj)
            self.model.obj = pyo.Objective(rule=lambda x: self._obj_rule(x,alpha), sense=pyo.maximize)
        else:
            raise AttributeError(" ".join(["Couldn't find any instance to update.",
                "Make sure to call fit() on a model before editing."]))

    def solve(self, solver="glpk", verbose=True):
        s = solver.lower()
        if s not in self._solvers:
            raise ValueError("Unsupported solver requested. Available solvers are: {}".format(self._solvers))
        
        sol = pyo.SolverFactory(s)
        if self._timeout is not None:
            # If specified, set solver timeout.
            #
            # Default timeout is one minute regardless of the solver being invoked.
            # Unsetting the time limit is highly discouraged, given the 
            # computational burden of exact algorithms (like this)
            # even for medium-sized instances of the allignment problem.
            if s == "cplex":
                sol.options["timelimit"] = self._timeout
            elif s == "glpk":
                sol.options["tmlim"] = self._timeout
            elif s == "gurobi":
                sol.options["TimeLimit"] = self._timeout

        print("Solving instance...", file=sys.stderr)
        
        self.solved_model = sol.solve(self.model, tee=verbose)
        self._sol_json = self.solved_model.json_repn()
        #print(self._sol_json)
        # Check solver status
        status = self._sol_json["Solver"][0]["Status"]
        if status == "ok":
            data = np.array([self.model.x[i,k]() for i in self.model.I for k in self.model.K])
            self.opt = csr_matrix(data.reshape((self._n1, self._n2)))
            self.opt_val = self.model.obj()
        else:
            raise RuntimeError("Solver status:", status)

    def node_similarity(self):
        bc = sum(self.S[(i,k)] * self.model.x[i,k]() for k in self.model.K for i in self.model.I)
 
        return bc / self._n1
    
    def edge_correctness(self):
        tc = sum(self.model.y[i,k]() for k in self.model.K for i in self.model.I)
        
        return tc / self.A.sum()

    def get_wtime(self):
        return self._sol_json["Solver"][0]["Time"]
