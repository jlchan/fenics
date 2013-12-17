# confusion for DPG with primal formulation

from dolfin import *

zero = Constant(0.0)
beta = Constant((1.0,.5))
eps = 1.0e-1

class Inflow(Expression):
	def eval(self, values, x):
		values[0] = 0.0
		if abs(x[0]) < 1E-14 or abs(x[1]) < 1E-14:
			values[0] = 1.0
inflowIndicator = Inflow()
outflowIndicator = Expression('1.0')-inflowIndicator
def inflow(x):
	return abs(x[0]) < 1E-14 or abs(x[1]) < 1E-14
def outflow(x):
	return abs(x[0]-1) < 1E-14 or abs(x[1]-1) < 1E-14
	

degree = 1             # trial space degree p
N =  16
mesh = UnitSquareMesh(N,N)

E   = FunctionSpace(mesh,"CG", degree+2)          # error estimator e
CG  = FunctionSpace(mesh,"CG", degree+1)          # primal variable u
Q   = FunctionSpace(mesh,"BDM",degree,restriction="facet") # interfacial flux q                                                 
X = MixedFunctionSpace( [E,CG,Q] )
(e,u,f) = TrialFunctions(X)
(v,du,df) = TestFunctions(X)

n = FacetNormal(mesh)                             # normal vectors

# inner product
ip = inner(dot(beta,grad(e)),dot(beta,grad(v)))*dx + eps*inner(grad(e),grad(v))*dx + inner(e,v)*dx 

def b(w,t,q):
	return inner(grad(w),beta*q + eps*grad(q)) * dx - inner(dot(t('+'),n('+')),(q('+')-q('-'))) * dS - inner(inflowIndicator*dot(t,n),v)*ds 

a = b(u,f,v) + b(du,df,e) + ip 

F = Expression("1.0")
L = F*v*dx                                       # rhs

# Dirichlet boundary condition on al boundary
bc1 = DirichletBC(X.sub(1), zero, DomainBoundary())
bcs = bc1
bc2 = DirichletBC(X.sub(0), zero, outflow)
#bcs.append(bc2)

# solve
uSol = Function(X) 
solve( a==L, uSol, bcs)
(e,u,q) = uSol.split()

fineSpace = FunctionSpace(refine(mesh), "Lagrange", 1)
uF = interpolate(u, fineSpace)
eF = interpolate(e, fineSpace)
#plot(uF)
plot(u)
plot(mesh)
interactive()

file = File('u.pvd')
file << uF
file = File('mesh.pvd')
file << mesh
file = File('e.pvd')
file << eF
