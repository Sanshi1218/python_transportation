'''__author__ = 'xulei'''




from __future__ import division
from pyomo.environ import *
from pyomo.opt.base import SolverFactory
from pyomo.opt.parallel import SolverManagerFactory
from pyomo.opt.parallel.manager import solve_all_instances

from numpy import *


from pathopt import *
import pathopt as pmod
import time
import setting as set
import copy

f = open('path_output.txt','w')





start_time = time.time()

def run_path(data):
    #call solver
    opt = SolverFactory('cplex')
    solver_manager = SolverManagerFactory('pyro')
    #report_timing = True to see time taken for each model component
    inst = pmod.model.create_instance(data, report_timing = True)
    #fix variable I to be 0
    for (h,v,t) in inst.aux:
        inst.I[h,v,t].fix(0)
    #tee=True to trace the cplex process
    result = opt.solve(inst, tee=True, warmstart=True)
    inst.solutions.load_from(result)
    '''
    for (s,e) in inst.OD:
         for p in inst.P:
             for j in inst.N:
                for t in inst.TS:
                    if inst.q[s,e,p,j,t]() > 0:
                        print 'queue', s,e,p,j,t,inst.q[s,e,p,j,t]()
    '''
    for (s,e) in inst.OD:
        for p in inst.P:
            for j in inst.N:
                for t in inst.TS:
                    if inst.f[s,e,p,j,t]() > 0:
                        print >> f, 'outflow', s,e,p,j,t,inst.f[s,e,p,j,t]()
    '''
    for (s,e) in inst.OD:
        for p in inst.P:
            for j in inst.N:
                for t in inst.TS:
                    for v in inst.VS:
                        if inst.b[s,e,p,v,t,j]() > 0:
                            print 'board', s,e,p,v,t,j,inst.b[s,e,p,v,t,j]()
    '''
    print >>f, 'objective value', inst.obj()
    #result.write()
    inst.display()
    exec_time = time.time() - start_time
    print "execution time", exec_time, 's'




run_path('data1.dat')



