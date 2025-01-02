from sympy import *

X, h = symbols('X h')
k, q = symbols('k q')

x = IndexedBase('x')
dx = IndexedBase('dx')
y = IndexedBase('y')


# Below the polynomial order is called N
Nmax = 5

##################################################################
# we use Lagrange characteristic polynomials l^N_q(X) of order N
##################################################################

# c = 1/denominator in l^N_q(X)
c = []
for N in range(Nmax+1):
    dd = []
    for i in range(N+1):
        cc = dx[q]/( Product((dx[k]), (k, 0, N)).doit() )
        cc = cc.replace(dx[q], dx[i])
        for r in range(N+1):
            cc = cc.replace(dx[r], (x[i]-x[r]))
        dd.append(cc)
    c.append(dd)

# numerator in l^N_q(X)
loc = []
for N in range(Nmax+1):
    ll = []
    for i in range(N+1):
        cc = ( Product((dx[k]), (k, 0, N)).doit() )/dx[q]
        cc = cc.replace(dx[q], dx[i])
        for r in range(N+1):
            cc = cc.replace(dx[r], (X-x[r]))
        ll.append(cc)
    loc.append(ll)

# derivs of loc
dloc = []
for N in range(Nmax+1):
    dl = []
    for i in range(N+1):
        if N > 0:
            dd = diff(loc[N][i], X, N-1)
        else:
            dd = loc[N][i]
        dl.append(dd)
    dloc.append(dl)

# derivs of loc on uniform grid: x[j] = j*h
dlocu = []
for N in range(Nmax+1):
    dl = []
    for i in range(N+1):
        dd = dloc[N][i]
        for j in range(N+1):
            dd = dd.replace(x[j], j*h)
        dl.append(dd)
    dlocu.append(dl)

# coeffs on uniform grid: x[j] = j*h
cu = []
for N in range(Nmax+1):
    dl = []
    for i in range(N+1):
        dd = c[N][i]
        for j in range(N+1):
            dd = dd.replace(x[j], j*h)
        dl.append(dd)
    cu.append(dl)

# derivs of l uniform grid: x[j] = j*h
dlu = []
for N in range(Nmax+1):
    dl = []
    for i in range(N+1):
        dd = dlocu[N][i] * cu[N][i]
        dd = dd.simplify()
        dl.append(dd)
    dlu.append(dl)

# polynom for derivs on uniform grid
dpu = []
for N in range(Nmax+1):
    pp = 0
    for i in range(N+1):
        pp += y[i] * dlu[N][i]
        pp = pp.simplify()
    dpu.append(pp)

# stencils on uniform grid
dpui = []
for N in range(Nmax+1):
    pp = []
    for j in range(N+1):
        dpu_Xj = dpu[N].replace(X, j*h)
        dpu_Xj = dpu_Xj.simplify()
        pp.append(dpu_Xj)
    dpui.append(pp)

# get increase factor of deriv-stencils for highest freq noise
# i.e. what we get for y[k] = +1,-1,+1,-1,+1,-1,...
fac = []
for N in range(Nmax+1):
    pp = []
    for j in range(N+1):
        dp = dpui[N][j]
        s = 1
        for n in range(N+1):
            dp = dp.replace(y[n], s)
            s = -s
        pp.append(dp)
    pm = pp[int(N/2)]
    for n in range(N+1):
        #pp[n] = pm/pp[n]
        #pp[n] = pp[n].evalf(17)
        pp[n] = abs(pp[n]/pm)
    fac.append(pp)

#############################################################
# Print stencils and factors
#############################################################
def get_stencil(dpui, N,i):
    st = dpui[N][i]
    stenc = []
    for m in range(N+1):
        sm = st
        for n in range(N+1):
            sm = sm.replace(y[n], KroneckerDelta(m,n))
        sm = sm.replace(h, 1)
        sm = sm.simplify()
        stenc.append(sm)
    return stenc


# print stencils
for N in range(Nmax+1):
    if N%2 == 0: continue
    print()
    print('/* One-sided ', N-1, 'th order derivative stencil */', sep='')
    print('double bsw',N-1,'[][',N+1,'] = {' , sep='')
    for j in range(N+1):
        stenc = get_stencil(dpui, N,j)
        print('  {' , sep='', end=' ')
        for m in range(N+1):
            print(stenc[m], sep='', end='')
            if m<N:
                print(',', sep='', end=' ')
        if j<N:
            print(' },',sep='')
        else:
            print(' }')
    print('};')
    print('/* factor by which we should divide deriv for diss near boundary */')
    print('double bf',N-1,'[] = { ' , sep='', end='')
    for j in range(N+1):
        print(fac[N][j], sep='', end='')
        if j<N:
            print(', ',sep='', end='')
        else:
            print(' };')
