from dolfin import *
import numpy as np

# Create mesh and define function space
N=32;
p = 2;
mesh = UnitSquareMesh(N, N)

V = FunctionSpace(mesh, "Lagrange", p)
TAU = VectorFunctionSpace(mesh, "CG", p)
E = V*TAU

class BC(Expression):
	def eval(self, values, x):
		tol = 1e-14
		values[0] = 0.0
		if abs(x[0]) < tol:
			values[0] = 1.0
bc = BC()
eps = .01
print "eps = ", eps
beta = Expression(('1.0','0.0'))

# Define variational problem
(v,tau) = TrialFunctions(E)
(dv,dtau) = TestFunctions(E)

u = Expression("1.0")
x = V.cell().x
n = FacetNormal(mesh)
#u = conditional(le( (x[0]-.5)**2 + (x[1]-.5)**2 ,  0.01), 1.0, 0.0) # discontinuous QOI data - ball of .1 around .5,.5
sigma = Expression(('0.0','0.0'))
uhat = 1.0
fhat = Expression("0.0")

#a = inner(div(eps*tau)-dot(beta,grad(v)),div(eps*dtau)-dot(beta,grad(dv)))*dx + inner(tau-grad(v),dtau-grad(dv))*dx + inner(v,dv)*dx
a = inner(div(eps*tau)-dot(beta,grad(v)),div(eps*dtau)-dot(beta,grad(dv)))*dx + inner(tau+grad(v),dtau+grad(dv))*dx + inner(v,dv)*dx
L = inner(u,eps*div(dtau)-dot(beta,grad(dv)))*dx + inner(sigma,dtau+grad(dv))*dx + bc*fhat*dv*ds + bc*uhat*dot(dtau,n)*ds 

# Compute solution
u = Function(E)
solve(a == L, u)
(v,tau) = u.split()
taux,tauy = tau.split()
# Plot solution
plot(v)
plot(taux)
plot(tauy)
interactive()