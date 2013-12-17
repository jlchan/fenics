from dolfin import *
from numpy import array, argsort, amax

import random


# Create mesh and define function space
N1=8
N2=8
h = min(1.0/N1,1.0/N2);
mesh = UnitSquareMesh(N1,N2)
boundary = BoundaryMesh(mesh,'exterior')

# Define Dirichlet boundary (x = 0 or x = 1)
# define spaces
def on_boundary(x):
	tol = 1E-14;
	return abs(x[0]-1)<tol or abs(x[1]-1)<tol or abs(x[0])<tol or abs(x[1])<tol

randomALE = False;
if (randomALE):
	C = .5*h
	for x in mesh.coordinates():
		if (on_boundary(x)):
			print("")
		else:
			d = sqrt((x[0]-.5)**2+(x[1]-.5)**2);
			if (d>h/2):
				x[0] += C*(random.random()-.5)#C*d;
				x[1] += C*(random.random()-.5)#C*d;

# Move mesh
mesh.move(boundary)
numRefs = 5
for refIndex in xrange(numRefs):
	V = FunctionSpace(mesh, "Lagrange", 3)
	# Define boundary condition
	u0 = Constant(0.0)
	bc = DirichletBC(V, u0, on_boundary)
	n = FacetNormal(mesh)	
	h = CellSize(mesh)
	h_avg = (h('+') + h('-'))/2
	
	# Define variational problem
	u = TrialFunction(V)
	v = TestFunction(V)
	f = Expression("exp(1.0/((x[0]-.5)*(x[0]-.5) + (x[1]-.5)*(x[1]-.5) + .1))")
	a = inner(grad(u), grad(v))*dx 
	L = f*v*dx

	# Compute solution
	u = Function(V)
	solve(a == L, u, bc)
	
		# evaluate element error indicator
	DG0 = FunctionSpace(mesh, "DG", 0) # element indicator function
	w = TestFunction(DG0)
	vf = interpolate(f,V)
	#M = w*h**2*f**2*dx 
	M = h**2*w*(-div(grad(u))-f)**2*dx + avg(f)*avg(w)*avg(h)*jump(grad(u),n)**2*dS
	
	cell_energy = assemble(M)

	# define adaptive strategy
	cell_markers = MeshFunction("bool", mesh, mesh.topology().dim())
	cell_markers.set_all(False)
	
	factor = .15
	# greedy refinement scheme
	maxE = amax(cell_energy.array())
	for c in cells(mesh):
		cell_markers[c] = cell_energy[c.index()] > factor*maxE
	
	mesh = refine(mesh,cell_markers)

mesh = refine(mesh)
fineMesh = refine(mesh)
fineSpace = FunctionSpace(fineMesh, "Lagrange", 1)
uF = interpolate(u,fineSpace)
# Save solution in VTK format
file = File("poisson.pvd")
file << u

# Plot solution
plot(uF)
plot(mesh)
interactive()
