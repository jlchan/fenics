""" ----------------------------------------------------------

This is a FEniCS implementation of a DPG method for the
Dirichlet problem

    Delta u = f    on UnitSquare
          u = g    on boundary,

and is part of notes for graduate lectures  introducing DPG
methods.

[Disclaimer: This file worked as of May 2013 with FEniCS version
1.2.0. It may not work in past or future versions!]

----------------------------------------Jay Gopalakrishnan  """


degree=1             # trial space degree p
m =  8               # m x m square mesh divided into triangles

from dolfin import *

# Expressions for u and f (for error and source computation)
U = Expression("x[0]*(1-x[0])*x[1]*(1-x[1])")
F = Expression("2*x[1]*(1-x[1])+2*x[0]*(1-x[0])")
## You can set other solutions by solely changing U and F above, e.g.:
# U = Expression("exp(pow(x[0],2)+pow(x[1],2))*sin(x[0]+x[1])")
# F = Expression("-2*exp(pow(x[0],2)+pow(x[1],2))*(sin(x[0]+x[1])*(1+2*(pow(x[0],2)+pow(x[1],2))) + 2*cos(x[0]+x[1])*(x[0]+x[1]))")

msh = UnitSquareMesh(m,m)
E   = FunctionSpace(msh,"DG", degree+2)          # error estimator e
CG  = FunctionSpace(msh,"DG", degree)          # primal variable u
Q   = FunctionSpace(msh,"RT",degree,restriction="facet")
                                                 # interfacial flux q
X = MixedFunctionSpace( [E,CG,Q] )
(e,u,q) = TrialFunctions(X)
(y,z,r) = TestFunctions(X)
n = FacetNormal(msh)                             # normal vectors

# the Y-inner product
yip = dot(grad(e),grad(y))*dx + e*y*dx

# b( (u,q), y ) = (grad u, grad y) - <q.n, y>
b1  = dot( grad(u),grad(y) ) * dx               \
    - dot(q('+'),n('+')) * (y('+')-y('-')) * dS \
    - dot(q,n)*y*ds - inner(dot(grad(u),n),y)*ds 

# b( (z,r), e)
b2 =  dot( grad(e),grad(z) ) * dx               \
    - (e('+')-e('-')) * dot(r('+'),n('+')) * dS \
    - e*dot(r,n)*ds - inner(dot(grad(z),n),e)*ds 

a = yip + b1 + b2                                # mixed form of DPG 
b = F*y*dx                                       # rhs

# Dirichlet boundary condition on al boundary
bc = DirichletBC(X.sub(1), U, DomainBoundary())
bc = DirichletBC(X.sub(0), U, DomainBoundary())

# solve
x = Function(X) 
solve( a==b, x, bcs=bc)
e,u,q = x.split(deepcopy=True)

# compute errors
fmsh=refine(refine(msh))
er = errornorm(U,u,norm_type='H1',degree_rise=3,mesh=fmsh)
print " Case of %d x %d mesh with degree %d:"% (m,m,degree)
print " H1-norm of  (u - uh)   = %15.12f" % er
print " Error estimator norm   = %15.12f" % \
    sqrt(assemble( (dot(grad(e),grad(e)) +e*e) * dx ))

# compute the H1 projection of U (a std Galerkin solution)
uu = TrialFunction(CG)
vv = TestFunction(CG)
aa = (dot(grad(uu),grad(vv) )+uu*vv) *dx
bb = (F+U)*vv*dx
bc = DirichletBC(CG, U, DomainBoundary())
pu = Function(CG)
solve( aa==bb, pu, bcs=bc)
erp = errornorm(U,pu,norm_type='H1',degree_rise=3,mesh=fmsh)
print " Error in H1 Projection = %15.12f" % erp
if abs(erp)>1.e-12:
    print " Error:Projection ratio = %15.12f" % (er/erp)
print " H1-norm of (uh - proj) = %15.12f" % \
    errornorm(u,pu,norm_type='H1',degree_rise=0)

# Issues:
#
#    If you put degree=2 then the H1 best approximation error comes
#    out LARGER than the error of a numerical solution in the same
#    space! This is mathematically impossible.  Is FEniCS losing some
#    double-precision digits in some routine?

