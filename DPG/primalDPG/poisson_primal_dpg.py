""" ----------------------------------------------------------

This is an implementation of the primal DPG method for the 
Poisson equation with -Del u = 1 on unit square and 
homogeneous boundary conditions.

----------------------------------------------------------  """

from dolfin import *

# Create mesh and define function space
degree = 1
delta_p = 2
N = 4
mesh = UnitSquareMesh(N,N)

u0 = Constant(0.0)

# define spaces
U = FunctionSpace(mesh, "CG", degree)					# trial functions
V = FunctionSpace(mesh, "DG", degree+delta_p)			# discontinuous test functions
F = FunctionSpace(mesh, "RT", degree, restriction="facet") 	# fluxes
E = MixedFunctionSpace([U,V,F])
(u,e,f) = TrialFunctions(E)
(du,v,df) = TestFunctions(E)
n = FacetNormal(mesh)

def ip(e,v):
	return inner(grad(e),grad(v)) + inner(e,v)

def b(u,v,f):
	return inner(grad(u),grad(v))*dx - inner(dot(f('+'),n('+')),jump(v)) * dS - inner(dot(f,n),v)*ds 

a = b(u,v,f) + b(du,e,df) + ip(e,v)*dx

L = inner(Constant(1.0),v)*dx

bcs = []
bcs.append(DirichletBC(E.sub(0), u0, DomainBoundary())) # boundary conditions on u

uSol = Function(E)
solve(a==L, uSol, bcs)
(u,e,f) = uSol.split()

plot(u)
plot(e)
plot(mesh)
interactive()