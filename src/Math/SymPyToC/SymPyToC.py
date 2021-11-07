# we need mostly sympy and very little from numpy
from sympy import *
import numpy as np

import itertools
import textwrap

###########################################################################
#
# assemble_all_EqnComponents :
#  * Eqs are input as a list of strings
#  * Eqs are then split into a list of 2 string lists for LHS and RHS
#    (by get_LHS_RHS_strlist)
#  * The 2 string lists are then converted into Indexed object lists
#    (Commands like :Decl or :Text are left as strings in the LHS)
#  * :Decl Commands first have just '.' as RHS
#  * All sums are expanded
#  * :Decl Commands get a RHS of this form
#     [ [A11,A12,...], [B11,B12,...], ...]
#     i.e. a list that contains lists of all components of the
#     declared tensors (as Indexed object)
#  * assemble_all_EqnComponents returns this for normal Eqs:
#     [ [ [A11,A12,...], [B11,B12,...], ...],   # LHS part
#       [ [rA11,...],    [rB11,...], ... ] ]    # RHS part
#    (and the above for the :Decl RHS parts)
#
# apply_symmetries_to_all_EqnComponents :
#  * apply symmetries to LHS and RHS => RHS simplify
#  * remove some unneeded LHS comps (and their corresponding RHS)
#  * apply symmetries to RHS parts of :Decl lines
#  * remove unneeded components in RHS parts of :Decl lines

###############################################################
# We mainly use sympy IndexedBase objects for our tensors here
# Example for other uses of sympy.IndexedBase :
#
# from sympy import *
# x = IndexedBase('x')
#j = Idx('j', (1, 3))
#k = Idx('k', (1, 3))
#n = Idx('n', (1, 3))
# f = 1/sqrt(Sum(x[k]**2, (k, 1, n)))
# print(f.diff(x[j]))


###########################################################################
# Functions set up tensor equations
###########################################################################

# convert tocompute string into [ [lhs0,lhs1, ...], [rhs0,rhs1, ...] ]
def get_LHS_RHS_strlist(tocompute):
    llist = []
    rlist = []
    for eqn in tocompute:
        if eqn[0] == ':':
            lhs = eqn
            rhs = '.'
        else:
            line = eqn.split('=', 1)
            lhs = line[0].strip()
            rhs = line[1] # line[1].strip()
        llist.append(lhs)
        rlist.append(rhs)
    return [llist,rlist]


# make IndexedBase objects from a string list
def make_IndexedObj_from_strlist(list):
    for vari in list:
        if vari[0] == ':':
            continue
        basestr = vari.split('[')[0]
        basestr = basestr.strip()
        #print(basestr,'|')
        cmd = basestr + ' = IndexedBase(\'' + basestr + '\')'
        print(' ', cmd)
        exec(cmd, globals())


#def recurse over contraction_structure to get final sum pieces
def sums_in_contraction_structure(ContrStruc):
    """
    This function returns the dict contres. contres contain the same keys
    dumm index keys as ContrStruc, whithout the keys that do not indicate
    summation (i.e. the ones that are just terms).
    The value corresponding to each key in contres, is itself a dict.
    This dict contains all the terms that needed to be summed as keys,
    and all the results of these partial sums as values. When we add
    these individual sums we get the result of the entire contraction:
    sum = 0
    for indkey in contres:
        termsdict = contres[indkey]
        for termkey in termsdict:
            sum += termsdict[termkey]

    To see what's in there we can print it:
    contres = sums_in_contraction_structure(ContrStruc)
    for indkey in contres:
        print('indkey =',indkey)
        termsdict = contres[indkey]
        for termkey in termsdict:
            print('termkey =',termkey)
    """
    #print('ContrStruc =', ContrStruc)

    # dict with results for each key
    contres = {} # empty dict

    # go over keys and sum or call sum_over_contraction_structure again
    for key in ContrStruc:
        #print('key =',key)

        # skip all keys that do not contain contraction indices
        if (key != None) and (type(key) != tuple):
            continue

        # assemble indicies we need to sum over
        indlist = []
        if type(key) == tuple:
            for ind in key:
                indlist.append(ind)

        contpartset = ContrStruc[key]
        contpartres = {} # empty dict

        #print('indlist =',indlist)
        #print('contpartset =',contpartset)

        # do all sub sums
        for contpart in contpartset:
            #print('contpart =',contpart)
            subres = contpart
            if contpart in ContrStruc:
                cslist = ContrStruc[contpart]
                sums_in_cs_list = []
                for cs in cslist:
                    sums_in_cs = sums_in_contraction_structure(cs)
                    sums_in_cs_list.append(sums_in_cs)
                #print('sums_in_cs_list =',sums_in_cs_list)

                # insert results in sums_in_cs_list into contpart -> subres
                for sums_in_cs in sums_in_cs_list:
                    for inds2 in sums_in_cs:
                        termsdict = sums_in_cs[inds2]
                        for termkey in termsdict:
                            res = termsdict[termkey]
                            #print('termkey =',termkey)
                            #print('res =',res)
                            subres = subres.subs(termkey, res)
                            #print('subres =',subres)
            # now use subres
            subp = subres
            #print('subp =',subp)
            for ind in indlist:
                # do sum over index ind and update subp
                subp = Sum(subp, (ind, ind.lower,ind.upper)).doit()

            # then add subp to the results dict
            contpartres[contpart] = subp

            contres[key] = contpartres

    return contres


#sum all terms is contraction_structure to get final sum
def sum_over_contraction_structure(ContrStruc):
    contres = sums_in_contraction_structure(ContrStruc)
    #print('sum_over_contraction_structure:')
    ##print('contres =',contres)
    #for indkey in contres:
    #    print('indkey =',indkey)
    #    termsdict = contres[indkey]
    #    for termkey in termsdict:
    #        print('termkey =',termkey)
    sum = 0
    for indkey in contres:
        termsdict = contres[indkey]
        for termkey in termsdict:
            sum += termsdict[termkey]
    return sum


# expand all sums in the term rhs
def expand_sums(rhs):
    tgi  = tensor.get_indices(rhs)
    tgcs = tensor.get_contraction_structure(rhs)
    sum = sum_over_contraction_structure(tgcs)
    return sum


# expand sums in all RHS
def expand_RHS_sums(eqs):
    LHS_list = []
    RHS_list = []
    for eq_i in range(len(eqs[0])):
        #print(i)
        ls = eqs[0][eq_i]
        rs = eqs[1][eq_i]
        #print(ls, rs)

        if ls[0] == ':':
            LHS_list.append(ls)
            RHS_list.append(rs)
            continue

        cmd = 'lhs = ' + ls
        #print(cmd)
        exec(cmd, globals())
        #print(lhs)

        cmd = 'rhs = ' + rs
        #print(cmd)
        exec(cmd, globals())
        #print(rhs)

        #print(lhs, '=', rhs)
        tgi_l  = tensor.get_indices(lhs)
        tgcs_l = tensor.get_contraction_structure(lhs)

        # check rhs is maybe just a constant number
        sympified_rhs = sympify(rhs)
        if sympified_rhs.is_Number:
            exprhs = rhs # do nothing
        else:
            tgi_r  = tensor.get_indices(rhs)
            tgcs_r = tensor.get_contraction_structure(rhs)
            #print(tgi_l, tgi_r)
            #print(tgcs_l, tgcs_r)
            exprhs = expand_sums(rhs) # do sums in RHS

        RHS_list.append(exprhs)
        LHS_list.append(lhs)
    return [LHS_list, RHS_list]


# return all possible index values of a tensor
def get_tuple_with_all_indexvals(T):
    if type(T) == str:
        return ()
    Tshape = T.shape
    if Tshape == None:
        return ()
    shap = tuple(Tshape) # maybe use:  tensor.get_indices(T)
    l1 = list(np.ndindex(shap))
    l2 = [ind.lower for ind in T.indices]
    a1 = np.array(l1, dtype=int)
    a2 = np.array(l2, dtype=int)
    a = list(a1 + a2)
    tup = tuple(a)
    return tup


# get expanded eqn components
def all_EqnComponents(expanded_eqs):
    # lists for all comps of LHS and RHS
    LHS_list = []
    RHS_list = []
    # loop over all Eqs
    for eqn_number in range(len(expanded_eqs[0])):
        # get all LHS components for this eqn_number
        lhs = expanded_eqs[0][eqn_number]
        lhs_allinds = get_tuple_with_all_indexvals(lhs)
        lhs_comps = [lhs.base[t] for t in lhs_allinds]
        if lhs_comps == []:
            lhs_comps = [lhs]
        # make subst rules for all indices
        subrules = []
        for indvals in lhs_allinds:
            sub = []
            for ind in range(len(indvals)):
                sub.append( (lhs.indices[ind], indvals[ind],) )
            subrules.append(sub)
        # find all RHS components for this eqn_number
        rhs = expanded_eqs[1][eqn_number]
        #print('lhs =',lhs)
        #print('rhs =',rhs)
        #print('subrules =',subrules)
        # make sure we can do subs on RHS
        if type(lhs) != str:
            sympified_rhs = sympify(rhs)
            if sympified_rhs.is_Number:
                # use same RHS for all components of LHS
                rhs_comps = [rhs for t in lhs_allinds]
            else:
                # get all RHS components for this eqn_number
                rhs_comps = [rhs.subs(s) for s in subrules]
        else:
            rhs_comps = [rhs]
        if rhs_comps == []:
            rhs_comps = [rhs]
        #print('rhs_comps =',rhs_comps)

        LHS_list.append(lhs_comps)
        RHS_list.append(rhs_comps)
    return [LHS_list, RHS_list]


# return all components of all Eqns as sympy objects
def assemble_all_EqnComponents(tocompute):
    #global AUTOVARS
    print('SymPyToC.py:')
    # print tensor Eqs
    #print('==========')
    #print('Equations:')
    #print('==========')
    #for s in tocompute:
    #    print(s)
    # make list of LHSs and RHSs of Eqns
    eqs = get_LHS_RHS_strlist(tocompute)
    print('Processing declarations')
    Declvars, LHSvars = Declvars_LHSvars_FromEqs(eqs)
    # introduce sympy vars
    make_IndexedObj_from_strlist(Declvars)
    make_IndexedObj_from_strlist(LHSvars)
    # construct the AUTOVARS
    print('Finding AUTOVARS')
    AUTOVARS = get_AUTOVARS(Declvars, LHSvars)
    # expand RHS of Decl parts
    print('Adding declarations')
    add_AllComponentsToDeclRHS(eqs, AUTOVARS)
    # expand sums on all RHS
    print('Expanding summations')
    expeqs = expand_RHS_sums(eqs)
    # get all components of all eqns
    print('Getting components')
    all = all_EqnComponents(expeqs)
    return all, AUTOVARS


# make AUTOVARS list, i.e. the vars that are in the LHS but not in Declvars
def get_AUTOVARS(Declvars, LHSvars):
    #print('Declvars =',Declvars)
    #print('LHSvars =',LHSvars)
    AUTOVARS1 = []
    # make list of Indexed objects from string list Declvars
    DeclTens = []
    for decv in Declvars:
        DeclTens.append(eval(decv))
    # loop over LHSs
    for lhs in LHSvars:
        Tl = eval(lhs)

        # check if Tl is in DeclTens
        lhs_is_in_DeclTens = False
        for decT in DeclTens:
            typ_decT = type(decT)
            typ_Tl   = type(Tl)
            if typ_decT != typ_Tl:
                continue
            if typ_Tl == tensor.indexed.Indexed:
                if Tl.base == decT.base:
                    if Tl.shape == decT.shape:
                        lhs_is_in_DeclTens = True
                        break
            else:
                if Tl == decT:
                    lhs_is_in_DeclTens = True
                    break

        if lhs_is_in_DeclTens == False:
            AUTOVARS1.append(lhs)

    # remove duplicates from AUTOVARS1
    AUTOVARS = []
    for av in AUTOVARS1:
        if av not in AUTOVARS:
            AUTOVARS.append(av)

    return AUTOVARS


###########################################################################
# Functions to deal with symmetries
###########################################################################

# count how many ordering defects a tuple has
def OrderingDefects(tup):
    ndefects = 0
    tprev = tup[0]
    for t in tup:
        #print('t =',t)
        if t < tprev:
            ndefects += 1
            #print('defect:', ndefects)
        tprev = t
    return ndefects

# find where lists like [i,j,k] and [i,k,j] differ
def findlistdiff(indices, sym):
    diff = []
    indnum = len(indices)
    for ind in range(indnum):
        if indices[ind] != sym[ind]:
            for ind2 in range(indnum):
                if sym[ind2] == indices[ind]: break
            #print(indices[ind], sym[ind], ind, ind2)
            # now we have ind and ind2
            diff.append( (ind, ind2) )
    return diff

# use to check if diff returned by findlistdiff in non-empty
def exit_if_list_empty(li, message):
    #print('li =',li)
    if li == []:
        raise SystemExit(message)

# return swapped indexlist
def swapped_indices(indexlist, diff):
    swpd = list(indexlist).copy()
    for tup in diff:
        swpd[tup[1]] = indexlist[tup[0]]
    return swpd

# return tensor comp with indices swaped if this improves ordering
def SwapIndicesIfOrderIsImproved(Tcomp, diff, symkind):
    #print(':', Tcomp, diff, symkind)
    indices = Tcomp.indices
    swpdindices = swapped_indices(indices, diff)
    #print(diff, indices, swpdindices)
    indsindiff = [ indices[t[0]] for t in diff]
    swapindiff = [ swpdindices[t[0]] for t in diff]
    #if Tcomp == T[1,2,2,1,3]:
    #    print('*', Tcomp, diff, symkind,':')
    #    print('  ', indices,'->',swpdindices,':',
    #          indsindiff,OrderingDefects(indsindiff),
    #          swapindiff,OrderingDefects(swapindiff))
    #    print('  ', indices,'->',swpdindices,':',
    #          OrderingDefects(indices),'->',
    #          OrderingDefects(swpdindices))
    # check what kind of symmetry we are dealing with:
    if(symkind == '+'):
        if(OrderingDefects(swapindiff) < OrderingDefects(indsindiff)):
            return 1, Tcomp.base[swpdindices]
    elif(symkind == '-'):
        if(OrderingDefects(swapindiff) <= OrderingDefects(indsindiff)):
            return -1, Tcomp.base[swpdindices]
    else:
            raise SystemExit('unknown symmetry kind')
    # default is to do nothing with Tcomp:
    return 1, Tcomp

# make substitution rules from symmetries
#TI is e.g. T[i,j,k,l,m]
def make_subsrules_from_symmetries(symmetries):
    subsruledict = {}
    for TI in symmetries:
        symsdict = symmetries[TI]
        # list with all symmetries
        slist = []
        for symkind in symsdict:
            for sym in symsdict[symkind]:
                #print(symkind, sym)
                slist.append( [symkind, sym] )
        # list with all permutations of all symmetries
        sym_perms = list(itertools.permutations(slist))
        # make all possible TI index combinations
        tup = get_tuple_with_all_indexvals(TI)
        allT = [TI.base[t] for t in tup] # here T is tensor head
        simpT = []                       # what we simplify it to
        for ci in range(len(allT)):
            comp = allT[ci]
            simpT.append([])
            #print(comp)
            #print(slist, sym_perms)
            comp_new = comp
            sign_new = 1
            # iterate as many times as there are permutations
            for iter in range(len(sym_perms)):
                #print('iter =', iter)
                for sp in sym_perms:
                    # try all minus symmetries and see if one gives zero
                    signfac = sign_new
                    cmp = comp_new
                    for s in sp:
                        symkind = s[0]
                        sym     = s[1]
                        if symkind == '-':
                            diff = findlistdiff(TI.indices, sym)
                            exit_if_list_empty(diff, 'error in symmetry of {T}'.format(T=TI))
                            cmp_old = cmp
                            sign, cmp = SwapIndicesIfOrderIsImproved(cmp,
                                                             diff, symkind)
                            # check if we have swapped on an antisymm pair
                            if sign==-1 and cmp == cmp_old:
                                signfac = 0
                                break
                    # try all syms in one sym_perm
                    sign_new = signfac
                    cmp = comp_new
                    for s in sp:
                        #print(sym_perms, s)
                        symkind = s[0]
                        sym     = s[1]
                        diff = findlistdiff(TI.indices, sym)
                        exit_if_list_empty(diff, 'error in symmetry of {T}'.format(T=TI))
                        # make sure this writes back into simpT:
                        cmp_old = cmp
                        sign, cmp = SwapIndicesIfOrderIsImproved(cmp,
                                                             diff, symkind)
                        # check if we have swapped on an antisymm pair
                        if sign==-1 and cmp == cmp_old:
                          signfac = 0
                        signfac *= sign
                    # if this symmetry improves or keeps order quality, we save it
                    if OrderingDefects(cmp.indices) <= OrderingDefects(comp_new.indices):
                        sign_new = signfac
                        comp_new = cmp
                    # if we already know it's zero save that
                    if signfac == 0:
                        sign_new = signfac
                        comp_new = cmp
            # we have tried all symmetry perms several times,
            # now save best result
            simpT[ci] = [sign_new, comp_new]
        # in subsrules we keep only those allT, simpT entries that differ
        subsrules = []
        for ind in range(len(allT)):
            if allT[ind] != simpT[ind][0] * simpT[ind][1]:
                subsrules.append( [allT[ind], simpT[ind]] )
        if subsrules != []:
            subsruledict[TI] = subsrules
    return subsruledict


###########################################################################
# Functions to simplify using the rules we derived from symmetries
###########################################################################

# use a subsruledict to simplify expr
def apply_subsrulesdict(subsruledict, expr):
    # do nothing for strings
    if type(expr) == str:
        return expr

    # recurse for lists
    if type(expr) == list:
        exprlist = expr
        for expr_i in range(len(exprlist)):
            ex_new = apply_subsrulesdict(subsruledict, exprlist[expr_i])
            exprlist[expr_i] = ex_new
        return exprlist

    # do nothing for numbers
    sympified_expr = sympify(expr)
    if sympified_expr.is_Number:
        return expr # do nothing

    # if we get here it should be a real math expression
    new_expr = expr
    for key in subsruledict:
        #if key.base != expr.base continue
        val = subsruledict[key]
        for sub in val:
            new_expr = new_expr.subs(sub[0], sub[1][0] * sub[1][1])
    return new_expr


# apply symmetries to all Eqn components
def apply_symmetries_to_all_EqnComponents(symmetries, allEqs, AUTOVARS):
    # remove Kdelta and LCeps3 first
    allEqs = simplify_all_EqnComponents(replace_Kdelta, allEqs)
    allEqs = simplify_all_EqnComponents(replace_LCeps3, allEqs)
    # make rules based on symmetries
    print('Making substitution rules from symmetries')
    subsruledict = make_subsrules_from_symmetries(symmetries)
    # use subsruledict to simplyfy LHS and RHS
    print('Applying symmetry substitution rules')
    for eq_i in range(len(allEqs[0])):
        for comp in range(len(allEqs[0][eq_i])):
            #print('L =', allEqs[0][eq_i])
            #print('R =', allEqs[1][eq_i])
            allEqs[0][eq_i][comp] = apply_subsrulesdict(subsruledict, allEqs[0][eq_i][comp])
            allEqs[1][eq_i][comp] = apply_subsrulesdict(subsruledict, allEqs[1][eq_i][comp])

    # make list of Eqs that we actually need
    print('Removing unneeded equations')
    simpLHS = []
    simpRHS = []
    for eq_i in range(len(allEqs[0])):
        # record which LHS side comps we may want to keep
        compstokeep = []
        for comp in range(len(allEqs[0][eq_i])):
            LHSstr = str(allEqs[0][eq_i][comp]);
            if(LHSstr[0] == '0' or LHSstr[0] == '-'):
                continue
            compstokeep.append(comp)
        # get rid of duplicates in LHS within compstokeep
        LHScomps = []
        RHScomps = []
        for nnn in compstokeep:
            if allEqs[0][eq_i][nnn] in LHScomps:
                continue
            LHScomps.append(allEqs[0][eq_i][nnn])
            RHScomps.append(allEqs[1][eq_i][nnn])
        simpLHS.append(list(LHScomps))
        simpRHS.append(list(RHScomps))

    # go over RHS of :Decl lines (for now they are all lists)
    print('Removing unneeded declarations')
    for eq_i in range(len(simpLHS)):
        LHScomp0 = simpLHS[eq_i][0]
        # if we have a :Decl command
        if type(LHScomp0) == str:
            if LHScomp0.startswith(':Decl:'):
                RHS = simpRHS[eq_i][0]
                # go over lists of lists in RHS
                newRHS = []
                for var in RHS:
                    pruned = remove_UnneededComps(var)
                    newRHS.append(pruned)
                simpRHS[eq_i][0] = newRHS
    # set result
    allEqsComps = [simpLHS, simpRHS]
    #print('===========')
    #print('Components:')
    #print('===========')
    #for eq_i in range(len(allEqsComps[0])):
    #    for comp in range(len(allEqsComps[0][eq_i])):
    #        print(allEqsComps[0][eq_i][comp], '=', allEqsComps[1][eq_i][comp])
    return allEqsComps


# Remove duplicates, zeros or comps with minus from a tensor component list.
# Such terms can occur after applying symmetries
def remove_UnneededComps(Tcomps):
    # record which comps we may want to keep
    keep1 = []
    for comp in Tcomps:
        compstr = str(comp);
        if(compstr[0] == '0' or compstr[0] == '-'):
            continue
        keep1.append(comp)
    # get rid of duplicates within keep1
    tokeep = []
    for comp in keep1:
        if comp in tokeep:
            continue
        tokeep.append(comp)
    return tokeep


###########################################################################
# Functions to simplify further
###########################################################################

# apply function simp to RHS of allEqs
def simplify_all_EqnComponents(simp, allEqs):
    print('Simplifying with', simp)
    allRHS = allEqs[1].copy()
    for eq_i in range(len(allRHS)):
        for compn in range(len(allRHS[eq_i])):
            rhs = allRHS[eq_i][compn]
            if type(rhs) != str and type(rhs) != list and simp != None:
                sympified_rhs = sympify(rhs)
                if not sympified_rhs.is_Number:
                    rhs = simp(rhs)
                allRHS[eq_i][compn] = rhs
    return [allEqs[0], allRHS]


# global Functions called POWER, POW2, POW3 to be used instead of Pow
POWER = symbols('POWER', cls=Function)
POW2  = symbols('POW2',  cls=Function)
POW3  = symbols('POW3',  cls=Function)

# replace the Pow function of sympy to get rid of all ** or pow in the
# output
def replace_Pow(expr):
    # get symbolic representation of expr, and replace Pow
    s = srepr(expr)
    s = s.replace('Pow', 'POWER')
    # get expr without Pow
    expr1 = eval(s)
    # use wild card symbol w to match POWER(w,2) and POWER(w,3)
    w = Wild('w')
    #print(expr1)
    expr2 = expr1.replace(POWER(w,2), POW2(w))
    #print(expr2)
    expr3 = expr2.replace(POWER(w,3), POW3(w))
    return expr3

# global IndexedBase objects for Kronecker delta and 3d Levi Civita symbol
Kdelta = IndexedBase('Kdelta')
LCeps3 = IndexedBase('LCeps3')

# replace Kdelta by KroneckerDelta
def replace_Kdelta(expr):
    w1 = Wild('w1')
    w2 = Wild('w2')
    expr2 = expr.replace(Kdelta[w1,w2], KroneckerDelta(w1,w2))
    return expr2

#replace LCeps3 by LeviCivita
def replace_LCeps3(expr):
    w1 = Wild('w1')
    w2 = Wild('w2')
    w3 = Wild('w3')
    expr2 = expr.replace(LCeps3[w1,w2,w3], LeviCivita(w1,w2,w3))
    return expr2


###########################################################################
# Functions translate Eqs into C or some other language
###########################################################################

# construct comp string from var indicies
def make_CompString(indices):
    # construct string that labels component
    comp = ''
    for ind in indices:
        comp += str(ind)
    return comp


# old version of make_DeclList
def make_DeclList__old(IndexedObjList, fstr):
    declist = []
    # init counters
    lastvarbase = None
    listindex    = 0
    varcompindex = 0
    varid    = -1
    # loop over all objects
    for obj in IndexedObjList:
        # set varbase, indices
        if type(obj) == tensor.indexed.Indexed:
            varbase = obj.base
            indices = obj.indices
        else:
            varbase = obj
            indices = []
        # reset some counters if we encounter a new variable
        if varbase != lastvarbase:
            lastvarbase = varbase
            varcompindex = 0
            varid += 1
        # construct string that labels component
        comp = make_CompString(indices);
        # string s for one declaration
        s = fstr.format(VAR=varbase, COMP=comp, LI=listindex, CI=varcompindex,
                        VARID=varid)
        declist.append(s)
        # increment some counters
        listindex    += 1
        varcompindex += 1
        lastvarbase = varbase
    return declist


# Make a list of variable declarations from IndexedObjList and format
# string fstr
# make_DeclList can be used like this:
"""
declist = make_DeclList([ [A[1,2], A[1,3], A[2,3]], [g[1,1], g[1,2]] ],
'double *{VAR}{COMP} = '\
'Vard( node, Vind(vlu,FOCCZ4->i_{VAR} + {CI} ); // ID={VARID} listpos={LI}')
for s in declist: print(s)
"""
def make_DeclList(List_of_IndexedObjLists, fstr):
    declist = []
    # init counters
    listindex    = 0
    varid        = 0
    # loop over all vars
    for var in List_of_IndexedObjLists:
        # loop over all component objects
        varcompindex = 0
        for obj in var:
            # set varbase, indices
            if type(obj) == tensor.indexed.Indexed:
                varbase = obj.base
                indices = obj.indices
            else:
                varbase = obj
                indices = []
            # construct string that labels component
            comp = make_CompString(indices);
            # string s for one declaration
            s = fstr.format(VAR=varbase, COMP=comp, LI=listindex, CI=varcompindex,
                            VARID=varid)
            declist.append(s)
            # increment some counters
            varcompindex += 1
            listindex    += 1
        # reset some counters if we encounter a new variable
        varid += 1
    return declist


# read out all Declared vars and LHS vars from eqs
def Declvars_LHSvars_FromEqs(eqs):
    Declvars = []
    LHSvars = []
    for eq_i in range(len(eqs[0])):
        # check if it has a command (that all start with ':')
        lhs = eqs[0][eq_i]
        if lhs[0] == ':':
            comstr = lhs
            if comstr.startswith(':Decl:'):
                sub = comstr.split(':')
                # list of vars is right after :Decl:
                varsstr = sub[2]
                # assume list of vars is separated by ;
                vars = varsstr.split(';')
                # remove spaces and append to Declvars
                for var in vars:
                    ss = var.strip()
                    if ss == 'AUTOVARS':
                        continue
                    if len(ss) > 0:
                        Declvars.append(ss)
        else:
            LHSvars.append(lhs)
    #print(Declvars, LHSvars)
    return Declvars, LHSvars


# put all components of vars after :Decl: string into a list on RHS
def add_AllComponentsToDeclRHS(eqs, AUTOVARS):
    # lists for all comps of LHS and RHS
    LHS_list = []
    RHS_list = []
    # loop over all Eqs
    for eqn_number in range(len(eqs[0])):
        lhs = eqs[0][eqn_number]
        if lhs[0] == ':':
            comstr = lhs
            if comstr.startswith(':Decl:'):
                sub = comstr.split(':')
                # list of vars is right after :Decl:
                varsstr = sub[2]
                # assume list of vars is separated by ;
                vars = varsstr.split(';')
                vars[0] = vars[0].strip()
                if vars[0] == 'AUTOVARS':
                    todeclare = AUTOVARS
                else:
                    todeclare = []
                    for var in vars:
                        ss = var.strip()
                        if len(ss) > 0:
                            todeclare.append(ss)
                # convert todeclare into Indexed Obj Tensors
                Tlist = []
                for var in todeclare:
                    Tensor = eval(var)
                    Tlist.append(Tensor)
                # clear RHS
                eqs[1][eqn_number] = []
                # next get all components
                for Tensor in Tlist:
                    allinds = get_tuple_with_all_indexvals(Tensor)
                    Tensor_comps = [Tensor.base[t] for t in allinds]
                    if Tensor_comps == []:
                        Tensor_comps = [Tensor]
                    # save Tensor_comps in RHS
                    eqs[1][eqn_number].append(Tensor_comps)
    return

# read Access fields of :Decl: to find the format for tensor output
def get_TensorOutputFormat(allEqs, AUTOVARS):
    outformat = []
    for eq_i in range(len(allEqs[0])):
        Lvar = allEqs[0][eq_i]
        for nnn in range(len(Lvar)):
            LHS = Lvar[nnn]
            # Commands
            if type(LHS) == str:
                comstr = LHS
                if LHS.startswith(':Decl'):
                    #print('Command',comstr)
                    ind = comstr.find('Access')
                    if ind == -1:
                        continue
                    text = comstr[ind:]
                    ind = comstr.find('=', ind)
                    text = comstr[ind+1:]
                    text = text.split(':')[0]
                    fstr = text.strip()
                    varlist = allEqs[1][eq_i][nnn]
                    for var in varlist:
                        for obj in var:
                            # set varbase, indices
                            if type(obj) == tensor.indexed.Indexed:
                                varbase = obj.base
                                indices = obj.indices
                            else:
                                varbase = obj
                                indices = []
                            # construct string that labels component
                            comp = make_CompString(indices);
                            # string s for this var component
                            s = fstr.format(VAR=varbase, COMP=comp)
                                            #, LI=listindex, CI=varcompinde
                                            #VARID=varid)
                            outformat.append([str(obj), s])
    return outformat


# go over final eqs and output them into a file
def write_Eqs(filename, allEqs, AUTOVARS):
    # first replace Pow in all expressions
    allEqs = simplify_all_EqnComponents(replace_Pow, allEqs)

    # construct outformat
    outformat = get_TensorOutputFormat(allEqs, AUTOVARS)
    #exit(9)

    # now do the output
    with open(filename, 'w') as f:
        for eq_i in range(len(allEqs[0])):
            Lvar = allEqs[0][eq_i]
            for comp in range(len(Lvar)):
                LHS = Lvar[comp]
                # Commands
                if type(LHS) == str:
                    comstr = LHS
                    if LHS.startswith(':Text'):
                        ind = comstr.find('=')
                        text = comstr[ind+1:]
                        # cut out leading '\n'
                        if text[0] == '\n':
                            text = text[1:]
                        f.write(text)
                    elif LHS.startswith(':Decl'):
                        #print('Command',comstr)
                        ind = comstr.find('DeclFunc')
                        if ind != -1:
                            text = comstr[ind:]
                            ind = comstr.find('=', ind)
                            text = comstr[ind+1:]
                            text = text.split(':')[0]
                            text = text.strip()
                            DeclFunc = eval(text)
                        else:
                            DeclFunc = make_DeclList
                        ind = comstr.find('Format')
                        if ind == -1:
                            continue
                        text = comstr[ind:]
                        ind = comstr.find('=', ind)
                        text = comstr[ind+1:]
                        if text[0] == '\n':
                            text = text[1:]
                        fstr = text
                        decs = DeclFunc(allEqs[1][eq_i][comp], fstr)
                        for line in decs:
                            f.write(line)
                # Eqs
                else:
                    RHS = allEqs[1][eq_i][comp]
                    RHSstr = str(RHS)
                    LHSstr = str(LHS)
                    for fmt in outformat:
                        RHSstr = RHSstr.replace(*fmt)
                        LHSstr = LHSstr.replace(*fmt)
                    RHSstr = textwrap.fill(RHSstr)
                    f.write(LHSstr)
                    f.write('\n=\n')
                    f.write(RHSstr)
                    f.write(';\n\n')
    return
