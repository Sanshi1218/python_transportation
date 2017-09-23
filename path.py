#'__author__ = 'xulei'


####with timetable adjustment

from pyomo.environ import *     #remember its pyomo.environ instead of just pyomo!!!


################################pure path model without timetable adjustment##############################

model = AbstractModel() #this is a class

#define network sets and params
model.n = Param(within=NonNegativeIntegers) #number of nodes, within= is used to specify domain
model.N = RangeSet(1,model.n)   #set of nodes from 1 to n

model.A = Set(within=model.N * model.N) #set of arcs
model.Ai = Set(within=model.N * model.N)  #interchange arcs

model.OD = Set(within=model.N * model.N)    #set of OD pairs
model.P = RangeSet(1,2) #number of feasibel paths for each OD
model.PP = Set(model.OD,model.P, within=model.N)
model.A1 = Set(model.OD,model.P)
model.N1 = Set(model.OD,model.P)     #interchange out
model.N2 = Set(model.OD,model.P)     #interchange in
model.NP = Set(model.P)     #set of nodes included in path p


#define timetable sets and params
model.T = Param(within=NonNegativeIntegers)
model.TS = RangeSet(1,model.T)
model.TSql = RangeSet(0,model.T)
model.V = Param(within=NonNegativeIntegers) #number of total trains
model.VS = RangeSet(1,model.V)
model.V1 = Param(within=NonNegativeIntegers)    #number of exising trains
model.VE = RangeSet(1,model.V1)
model.V2 = Param(within=NonNegativeIntegers)
model.VA = RangeSet(model.V2,model.V) #set of added trains in future
model.L = Param(within=NonNegativeIntegers)
model.LS = RangeSet(1,model.L)
model.LN = Set(model.LS, within=model.N)  #set of nodes for line h
model.C = Param(model.VS, default=900, mutable=True)
model.D = Param(model.OD,model.N,model.TS, default=0)
model.tr = Param(model.LS, model.N) #cum travel time
model.Y = Set(within=model.VE * model.TS * model.N)     #current timetable (v,t,j)


###TIMETABLE ADJUSTMENT SET
model.KS = Set(model.LS, within=model.TS)
model.KS1 = Set(model.LS, within=model.TS)
model.KSY = Set(model.LS, within=model.TS * model.N)



#####define indicator functions
def X0_rule(model,v,t,j):
    if (v,t,j) in model.Y:
        return 1
    return 0
model.X0 = Param(model.VE, model.TS, model.N, initialize=X0_rule)

#Xf use also X1
def X1_rule(model,s,e,j):
    if j==e:
        return 1
    return 0
model.X1 = Param(model.OD, model.N, initialize=X1_rule)

def X2_rule(model,s,e,j):
    if j==s:
        return 1
    return 0
model.X2 = Param(model.OD, model.N, initialize=X2_rule)

def X3_rule(model,i,j):
    if (i,j) in model.Ai:
        return 1
    return 0
model.X3 = Param(model.A, initialize=X3_rule)

def Xb_rule(model,s,e,p,j):
    if j==s or j in model.N2[s,e,p]:
        return 1
    return 0
model.Xb = Param(model.OD, model.P, model.N, initialize=Xb_rule)

def Xa_rule(model,s,e,p,j):
    if j==e or j in model.N1[s,e,p]:
        return 1
    return 0
model.Xa = Param(model.OD, model.P, model.N, initialize=Xa_rule)

'''
def Xd_rule(s,e,p,i,j):
    if (i,j) in model.A1[s,e,p]:
        return 1
    return 0
model.Xd = Param(model.OD, model.P, model.A, initialize=Xd_rule)
'''



#define vars
model.q = Var(model.OD, model.P, model.N,model.TSql, domain = NonNegativeReals)
model.f = Var(model.OD, model.P, model.N, model.TS, domain = NonNegativeReals)
model.d = Var(model.OD, model.P, model.A, model.TS, domain = NonNegativeReals)
model.b = Var(model.OD, model.P, model.VS,model.TS, model.N, domain = NonNegativeReals)
model.a = Var(model.OD, model.P, model.VS,model.TS, model.N, domain = NonNegativeReals)
model.l = Var(model.OD, model.P, model.VS,model.TSql, domain = NonNegativeReals)
model.w = Var(model.OD, model.P, model.N, model.TS, domain = NonNegativeReals)


'''
#define param over an indexed set
def set_init(model):
  """Initializes auxiliary set."""
  return [(w, ll, sw) for w in model.W for ll in model.LL[w] for sw in model.SW[w]]
model.AUX = cpr.Set(dimen=3, initialize=set_init)
model.Load = cpr.Var(model.AUX, within=cpr.Binary)
'''


def set_init(model):
    return [(h,v,k) for h in model.LS for v in model.VA for k in model.KS[h]]
model.aux = Set(dimen=3, initialize=set_init)
model.I = Var(model.aux, domain=Binary)



#objective and constraints
def obj_rule(model):
     return sum(model.q[s,e,p,j,t]
                for (s,e) in model.OD
                for p in model.P
                for j in model.N
                for t in model.TS)
model.obj = Objective(rule = obj_rule, sense = minimize)

#####Most constraints are specified using equality or inequality expressions that are created using a ####rule####, which is a Python function.
def con_1_rule(model,s,e,p,j):
     return model.q[s,e,p,j,0] == 0
model.con1 = Constraint(model.OD, model.P, model.N, rule=con_1_rule)

def con_2_rule(model,s,e,p,j,t):     ###############need to specify that the set of these indexed belong to?????##################
     return model.q[s,e,p,j,t] == model.q[s,e,p,j,t-1] + model.w[s,e,p,j,t] + sum(model.a[s,e,p,v,t,j] for v in model.VS) - sum(model.b[s,e,p,v,t,j]
              for v in model.VS) - sum(model.d[s,e,p,j,i,t] for (j,i) in model.A) + sum(model.d[s,e,p,i,j,t] for (i,j) in model.A) - model.f[s,e,p,j,t]
model.con2 = Constraint(model.OD, model.P, model.N, model.TS, rule=con_2_rule)

def con_3_rule(model,s,e,p,v):
     return model.l[s,e,p,v,0] == 0
model.con3 = Constraint(model.OD, model.P, model.VS, rule=con_3_rule)
#
def con_4_rule(model,s,e,p,v,t):
     return model.l[s,e,p,v,t] == model.l[s,e,p,v,t-1] - sum(model.a[s,e,p,v,t,j] for j in model.N) + sum(model.b[s,e,p,v,t,j] for j in model.N)
model.con4 = Constraint(model.OD, model.P, model.VS, model.TS, rule=con_4_rule)

def con_5_rule(model,v,t):
    return sum(model.l[s,e,p,v,t]
               for (s,e) in model.OD
               for p in model.P) <= model.C[v]
model.con5 = Constraint(model.VS, model.TS, rule=con_5_rule)

def con_6_rule(model,s,e,j,t):
    return sum(model.w[s,e,p,j,t]
               for p in model.P) == model.D[s,e,j,t]
model.con6 = Constraint(model.OD, model.N, model.TS, rule=con_6_rule)

def con_7_rule(model,s,e,p,j,t):
     return model.w[s,e,p,j,t] <= model.D[s,e,j,t]
model.con7 = Constraint(model.OD, model.P,model.N, model.TS,rule = con_7_rule)
#
def con_8_rule(model,s,e,p,j,t):
    return model.q[s,e,p,j,t] <= (1-model.Xa[s,e,p,j]) * sum(model.D[s,e,s,t1] for t1 in model.TS if t1<=t)
model.con8 = Constraint(model.OD, model.P, model.N, model.TS, rule = con_8_rule)

def con_9_rule(model, s,e,p,i,j,t):
     if (i,j) in model.A1[s,e,p]:
         return model.d[s,e,p,i,j,t] <= sum(model.D[s,e,s,t1] for t1 in model.TS if t1<=t)
     else:
         return model.d[s,e,p,i,j,t] == 0
model.con9 = Constraint(model.OD, model.P, model.A, model.TS,rule = con_9_rule)

def con_10_rule(model,s,e,p,j,t):
    return model.f[s,e,p,j,t] <= model.X1[s,e,j]*sum(model.D[s,e,s,t1] for t1 in model.TS if t1<=t)
model.con10 = Constraint(model.OD, model.P, model.N, model.TS, rule = con_10_rule)

def con_11_rule(model, s,e,p,v,t,j):
    return model.b[s,e,p,v,t,j] <= model.X0[v,t,j] * model.Xb[s,e,p,j] * sum(model.D[s,e,s,t1] for t1 in model.TS if
                                                                           t1<=t)
model.con11 = Constraint(model.OD, model.P, model.VE, model.TS, model.N, rule = con_11_rule)

def con_14_rule(model, s,e,p,v,t,j):
    return model.a[s,e,p,v,t,j] <= model.X0[v,t,j] * model.Xa[s,e,p,j] * sum(model.D[s,e,s,t1] for t1 in model.TS if
                                                                              t1<=t)
model.con14 = Constraint(model.OD, model.P, model.VE, model.TS, model.N, rule = con_14_rule)

def con_17_rule(model,s,e,p,v,t,j):
    return model.l[s,e,p,v,t] <= (1-model.Xa[s,e,p,j]) * model.C[v]
def set_17(model):
    return [(s,e,p,v,t,j) for (s,e) in model.OD for p in model.P for v in model.VE for t in model.TS for j in model.N
            if (v,t,j) in model.Y]
model.aux17 = Set(dimen=6, initialize=set_17)
model.con17 = Constraint(model.aux17, rule = con_17_rule)



#######################################model.KS1, model.LN, or what #####################
#######################################model.KS1, model.LN, or what #####################
def con_12_rule(model, s,e,p,h,v,t,j,t1):
    return model.b[s,e,p,v,t,j] <= model.I[h,v,t1] * model.Xb[s,e,p,j] * sum(model.D[s,e,s,t2] for t2 in model.TS if
                                                                              t2 <= t)
def set_12(model):
    return [(s,e,p,h,v,t,j,t1) for (s,e) in model.OD for p in model.P for h in model.LS for v in model.VA for t in
            model.KS1[h] for j in model.LN[h] for t1 in model.KS[h] if t-t1==model.tr[h,j]]
model.aux12 = Set(dimen=8, initialize=set_12)
model.con12 = Constraint(model.aux12, rule = con_12_rule)

def con_13_rule(model, s,e,p,h,v,t,j):
    return model.b[s,e,p,v,t,j] == 0
def set_13(model):
    return [(s,e,p,h,v,t,j) for (s,e) in model.OD for p in model.P for h in model.LS for v in model.VA for t in
            model.TS for j in model.LN[h] if (t,j) not in model.KSY[h]]
model.aux13 = Set(dimen=7, initialize=set_13)
model.con13 = Constraint(model.aux13,rule = con_13_rule)


def con_15_rule(model,s,e,p,h,v,t,j,t1):
    return model.a[s,e,p,v,t,j] <= model.I[h,v,t1] * model.Xa[s,e,p,j] * sum(model.D[s,e,s,t2] for t2 in model.TS if
                                                                              t2 <= t)
def set_15(model):
    return [(s,e,p,h,v,t,j,t1) for (s,e) in model.OD for p in model.P for h in model.LS for v in model.VA for t in
            model.KS1[h] for j in model.LN[h] for t1 in model.KS[h] if t-t1==model.tr[h,j]]
model.aux15 = Set(dimen=8, initialize=set_15)
model.con15 = Constraint(model.aux15, rule = con_15_rule)

def con_16_rule(model, s,e,p,h,v,t,j):
    return model.a[s,e,p,v,t,j]==0
def set_16(model):
    return [(s,e,p,h,v,t,j) for (s,e) in model.OD for p in model.P for h in model.LS for v in model.VA for t in
            model.TS for j in model.LN[h] if (t,j) not in model.KSY[h]]
model.aux16 = Set(dimen=7, initialize=set_16)
model.con16 = Constraint(model.aux16,rule = con_16_rule)


def con_18_rule(model, s,e,p,h,v,t,j,t1):
    return model.l[s,e,p,v,t] <= (1-model.I[h,v,t1]) * 900
def set_18(model):
    return [(s,e,p,h,v,t,j,t1) for (s,e) in model.OD for p in model.P for h in model.LS for v in model.VA for t in
            model.TS for j in model.LN[h] for t1 in model.KS[h] if t-t1==model.tr[h,j] and (j in model.N1[s,e,
                                                                                                          p] or j==e)]
model.aux18 = Set(dimen=8, initialize=set_18)
model.con18 = Constraint(model.aux18, rule = con_18_rule)


def con_19_rule(model, s,e,p,j,t):
    return model.q[s,e,p,j,t] == 0
def set_19(model):
    return [(s,e,p,j,t) for (s,e) in model.OD for p in model.P for j in model.LN for t in model.TSql if j not in
            model.PP[s,e,p]]
model.aux19 = Set(dimen=5, initialize=set_19)
model.con19 = Constraint(model.aux19, rule = con_19_rule)


#######timetable#####
def con_20_rule(model,h,t):
    for t in model.KS[h]:
        return sum(model.I[h,v,t]
               for v in model.VA) <= 1
def set_20(model):
    return [(h,t) for h in model.LS for t in model.KS[h]]
model.aux20 = Set(dimen=2, initialize=set_20)
model.con20 = Constraint(model.aux20, rule=con_20_rule)

def con_21_rule(model,v):
    return sum(model.I[h,v,t]
               for h in model.LS
               for t in model.KS[h]) == 1
model.con21 = Constraint(model.VA, rule=con_21_rule)



