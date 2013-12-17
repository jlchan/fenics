from dolfin import *
from numpy import array, argsort, amax
from math import ceil

# Create mesh and define function space
pU = 1
pV = 2
N = 8
#mesh = UnitSquare(N,N)
mesh = UnitIntervalMesh(N)

# define problem params
eps = 1e-3
print "eps = ", eps
#beta = Expression(('1.0','1.0'))
beta = Expression('1.0')


u0 = Expression('0.0')
def u0_boundary(x, on_boundary):
	return on_boundary

class Inflow(Expression):
	def eval(self, values, x):
		values[0] = 0.0
		if abs(x[0]) < 1E-14:
			values[0] = 1.0
inflow = Inflow()

def outflow(x):
	return abs(x[0]-1) < 1E-14

# define spaces
U = FunctionSpace(mesh, "Lagrange", pU)
V = FunctionSpace(mesh, "Lagrange", pV)
E = U*V
(u,e) = TrialFunctions(E)
(du,v) = TestFunctions(E)
n = FacetNormal(mesh)
h = CellSize(mesh)

bc1 = DirichletBC(E.sub(0), u0, u0_boundary) # boundary conditions on u
bc2 = DirichletBC(E.sub(1), u0, outflow) # error boundary condition
#bc2 = DirichletBC(E.sub(1), Expression('0.0'), u0_boundary) # error boundary condition
bcs = [bc1, bc2]

a = inner(u.dx(0),v)*dx 
a = a + eps*inner(grad(u),grad(v))*dx - eps*inner(inflow*dot(grad(u),n),v)*ds
a = a + inner(du.dx(0),e)*dx + eps*inner(grad(du),grad(e))*dx - inner(eps*inflow*dot(grad(du),n),e)*ds # orthogonality condition
ip = inner(e.dx(0),v.dx(0))*dx + eps*inner(grad(e),grad(v))*dx + inner(e,v)*dx #+ eps*(pV**2/h)*inner(inflow*e,v)*ds
a = a + ip # adding residual contribution

f = Expression('1.0')
L = inner(f,v)*dx

uSol = Function(E)
solve(a==L, uSol, bcs)
(u,e) = uSol.split()

file = File('u1D.pvd')
file << u
file = File('e1D.pvd')
file << e
