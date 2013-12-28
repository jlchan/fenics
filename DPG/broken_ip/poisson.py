# confusion for DPG with primal formulation

from dolfin import *
import helper_functions as help

p = int(help.parseArg('--p',1))
fp = int(help.parseArg('--fp',p))
numRefs = 4
init_N = 4
h_vec = []
conf_err = []
nconf_err = []
for ref in xrange(0,numRefs):
    N = init_N*2**ref
    mesh = UnitSquareMesh(N,N)
    n = FacetNormal(mesh)
    V  = FunctionSpace(mesh,"CG", p)
    
    ue = Expression("sin(x[0]*pi)*sin(x[1]*pi)",pi = 2.0*acos(0.0))
    forcing = Expression("2.0*pi*pi*sin(pi*x[0])*sin(pi*x[1])",pi=2.0*acos(0.0)) + ue
    
    # ================ standard solve ====================
    
    u,v = TrialFunction(V),TestFunction(V)
    a = inner(grad(u),grad(v))*dx + inner(u,v)*dx
    L = forcing*v*dx # rhs
    bc = DirichletBC(V, Constant(0.0), DomainBoundary())
    uh = Function(V) 
    solve( a==L, uh, bc)
    
    # ================ nonconforming solve ==================
    
    Vh  = FunctionSpace(mesh,"DG", p)
    F = FunctionSpace(mesh,"RT",fp,restriction="facet") # interfacial flux
    #F = FunctionSpace(mesh,"CG",p,restriction="facet") # interfacial flux   
    #F = FunctionSpace(mesh,"BDM",p-1,restriction="facet") # interfacial flux - p RT and p-1 BDM appear to be equivalent
    
    E = Vh*F
    (u,f),(v,df) = TrialFunctions(E),TestFunctions(E)
    ah = inner(grad(u),grad(v))*dx + inner(u,v)*dx
    ah += inner(dot(f('+'),n('+')),v('+')-v('-')) * dS + inner(dot(f,n),v)*ds
    ah += inner(dot(df('+'),n('+')),u('+')-u('-')) * dS + inner(dot(df,n),u)*ds
    L = forcing*v*dx # rhs

    uf = Function(E)
    bc = DirichletBC(E.sub(0), Constant(0.0), DomainBoundary(), "pointwise")
    solve( ah==L, uf, bc)
    (uh_nc,f) = uf.split()
    
    # ================= error ========================

    errSpace = FunctionSpace(refine(mesh),"DG",p+2)

    e = ue - interpolate(uh,errSpace)
    l_err = sqrt(assemble(e**2*dx))
    print "l2 err = ", l_err
    conf_err.append(l_err)
    

    e_nc = ue - interpolate(uh_nc,errSpace)
    l_nc_err = sqrt(assemble(e_nc**2*dx))
    print "l2 non-conforming err = ", l_nc_err
    nconf_err.append(l_nc_err)
    h_vec.append(1.0/N)
    

r_c = help.compute_rates(h_vec,conf_err)
r_nc = help.compute_rates(h_vec,nconf_err)

print "h = ", h_vec
print "conf_err = ", conf_err
print "conf_rates = ", r_c
print "nconf_err = ", nconf_err
print "nconf_rates = ", r_nc

# ================= plotting ============================== 

#plot(uh)
#plot(uh_nc)
#interactive()

