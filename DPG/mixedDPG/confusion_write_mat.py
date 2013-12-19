from dolfin import *
from numpy import array, argsort, amax
from math import ceil
import helper_functions as help

# Create mesh and define function space
pU = 1
pV = 2
N = 8
eps = 0;
mesh = UnitSquareMesh(N,N)

# define problem params
beta = Expression(('1.0','.70'))

u0 = Expression('0.0')
zero = Expression('0.0')

def inflow(x):
	return abs(x[0]) < 1E-14 or abs(x[1]) < 1E-14
def outflow(x):
	return abs(x[0]-1) < 1E-14 or abs(x[1]-1) < 1E-14
class Inflow(Expression):
	def eval(self, values, x):
		values[0] = 0.0
		if abs(x[0]) < 1E-14 or abs(x[1]) < 1E-14:
			values[0] = 1.0
inflowIndicator = Inflow()
outflowIndicator = 1-inflowIndicator;

# define spaces
U = FunctionSpace(mesh, "Lagrange", pU)
V = FunctionSpace(mesh, "Lagrange", pV)
E = U*V
(u,e) = TrialFunctions(E)
(du,v) = TestFunctions(E)
n = FacetNormal(mesh)

#bcs = [DirichletBC(E.sub(0), u0, inflow)] # boundary conditions on u
#bcs.append(DirichletBC(E.sub(1), zero, outflow)) # error boundary condition
#if eps>1e-8:
#	bcs.append(DirichletBC(E.sub(0), u0, outflow)) # boundary conditions on u

def ip(e,v):
	return inner(dot(beta,grad(e)),dot(beta,grad(v))) + inner(e,v)
#	return inner(dot(beta,grad(e)),dot(beta,grad(v))) + eps*inner(grad(e),grad(v)) + inner(e,v)

def b(u,v):
	return inner(-u,dot(beta,grad(v)))*dx + inner(outflowIndicator*dot(beta,n)*u,v)*ds 
#	return inner(grad(u),beta*v + (eps)*grad(v))*dx - inner(eps*inflowIndicator*dot(grad(u),n),v)*ds 

a = b(u,v) + b(du,e) + ip(e,v)*dx

f = Expression('1.0') 
L = inner(f,v)*dx
Q = Expression('x[0]*(1-x[0])*x[1]*(1-x[1])')
Lq = inner(Q,du)*dx;

# ==================================================

uSol = Function(E)
#solve(a==L, uSol, bcs)
solve(a==L, uSol)
(u,e) = uSol.split()
plot(u)
interactive()

A = assemble(a)
b = assemble(Lq)

#bcs[0].apply(A,b)
#bcs[1].apply(A,b)
#if eps>1e-8:
#	bcs[2].apply(A,b)

AA = A.array()
bb = b.array()

Udofs = E.sub(0).dofmap().collapse(mesh)[1].values()
Vdofs = E.sub(1).dofmap().collapse(mesh)[1].values()

uu = TrialFunction(U)
vv = TestFunction(U)
precondRiesz = inner(uu,vv)*dx #+ eps*inner(grad(uu),grad(vv))*dx
P = assemble(precondRiesz)
#trash = assemble(inner(Q,vv)*dx)
#bcU = [DirichletBC(U, u0, inflow), DirichletBC(U, u0, outflow)];
#bcU[0].apply(P, trash)
#bcU[1].apply(P, trash)
PP = P.array()

#print "UDofs = ", Udofs
#print "VDofs = ", Vdofs
#print "A = ", AA

dumpMat("readA.m",AA)
dump_d_vec("readb.m",bb)
dump_i_ec("readUi.m",Udofs)
dump_i_ec("readVi.m",Vdofs)
dumpMat("readP.m",PP)
