from dolfin import *

useBulkChasing = True
numRefs = 4

# Create mesh and define function space
pU = 1
pV = 3
N = 16
mesh = UnitSquareMesh(N,N)

# define problem params
eps = 1e-2
print "eps = ", eps
beta = Expression(('.5','1.0'))

u0 = Expression('0.0')
zero = Constant('0.0')
#u0 = Expression('sin(x[0]*pi*4.0)')
#u0 = Expression('(1-x[0])*(1-x[1])')
# restricts to leftInflow only
class Inflow(Expression):
	def eval(self, values, x):
		values[0] = 0.0
		if abs(x[0]) < 1E-14 or abs(x[1]) < 1E-14:
			values[0] = 1.0
inflowIndicator = Inflow()
outflowIndicator = Expression('1.0')-inflowIndicator

def u0_boundary(x, on_boundary):
	return on_boundary
def inflow(x):
	return abs(x[0]) < 1E-14 or abs(x[1]) < 1E-14
def outflow(x):
	return abs(x[0]-1) < 1E-14 or abs(x[1]-1) < 1E-14
class U0(Expression):
	def eval(self, values, x):
		tol = 1e-14
		values[0] = 0.0
		if abs(x[1]) < tol and x[0]<.5:
			values[0] = 1.0
		if abs(x[0]) < tol:
			values[0] = 1.0-x[1]
#u0 = U0()

#gaussian = Expression(".1*exp(-pow(x[0]-.5,2)*1e5)")

enrgy=[]
nDofs=[]
for refIndex in xrange(numRefs):
	# define spaces
	U = FunctionSpace(mesh, "Lagrange", pU)
#	Vb = FunctionSpace(mesh, "B", pV)
#	Vc = FunctionSpace(mesh, "Lagrange", pU)
#	V = Vc+Vb
	V = FunctionSpace(mesh, "Lagrange", pV)
	E = U*V
	(u,e) = TrialFunctions(E)
	(du,v) = TestFunctions(E)
	n = FacetNormal(mesh)
	h = CellSize(mesh)

	bcs = [DirichletBC(E.sub(0), u0, inflow)] # boundary conditions on u
	bcs.append(DirichletBC(E.sub(0), zero, outflow)) # boundary conditions on u
#	bcs.append(DirichletBC(E.sub(1), zero, inflow)) # error boundary condition
	bcs.append(DirichletBC(E.sub(1), zero, outflow)) # error boundary condition

	def ip(e,v):
		return inner(dot(beta,grad(e)),dot(beta,grad(v)))*dx + eps*inner(grad(e),grad(v))*dx + inner(eps*inflowIndicator*e,v)*ds
#		return inner(dot(beta,grad(e)),dot(beta,grad(v)))*dx + eps*inner(grad(e),grad(v))*dx + inner(e,v)*dx
	
	def b(u,v):
		return inner(grad(u),beta*v + (eps)*grad(v))*dx - inner(eps*inflowIndicator*dot(grad(u),n),v)*ds
#		return inner(grad(u),beta*v + (eps)*grad(v))*dx - inner(eps*inflowIndicator*dot(grad(u),n),v)*ds - inner(eps*outflowIndicator*dot(grad(v),n),u)*ds 

	a = b(u,v) + b(du,e) + ip(e,v)

	f = Expression('0.0')
	x = V.cell().x
	f = conditional(le( x[1]-2.0*x[0],  0.0), 1.0, 0.0) # discontinuous forcing data
	L = inner(f,v)*dx
	
	uSol = Function(E)
	solve(a==L, uSol, bcs)
	(u,e) = uSol.split()

	# evaluate element error indicator
	DG0 = FunctionSpace(mesh, "DG", 0) # element indicator function
	w = TestFunction(DG0)
	M = w*inner(dot(beta,grad(e)),dot(beta,grad(e)))*dx + w*eps*inner(grad(e),grad(e))*dx
	cell_energy = assemble(M)	

	# define adaptive strategy
	cell_markers = MeshFunction("bool", mesh, mesh.topology().dim())
	cell_markers.set_all(False)

	factor = .25
	if useBulkChasing:
		cutoff = sorted(cell_energy, reverse=True)[int(len(cell_energy)*factor)]
		for c in cells(mesh):
			cell_markers[c] = cell_energy[c.index()] > cutoff
	else:	
		# greedy refinement scheme
		maxE = amax(cell_energy.array())
		for c in cells(mesh):
			cell_markers[c] = cell_energy[c.index()] > factor*maxE
	
	mesh = refine(mesh,cell_markers)
#	mesh = refine(mesh)
	print "on refinement ", refIndex

#	plot(u)
#	plot(mesh)

	energy = inner(dot(beta,grad(e)),dot(beta,grad(e)))*dx + eps*inner(grad(e),grad(e))*dx
	totalE = assemble(energy)
	enrgy.append(totalE)
	nDofs.append(U.dofmap().global_dimension())
#	print "Total energy = ", totalE	

from math import log as ln  # (log is a dolfin name too - and logg :-)
print 'energy(',1,') = ', enrgy[0], '; ndofs(',1,')=',nDofs[0]
for i in range(1, len(enrgy)):
	print 'energy(',i+1,') = ', enrgy[i],'; ndofs(',i+1,')=',nDofs[i]

uSol = Function(E)
solve(a==L, uSol, bcs)
(u,e) = uSol.split()

fineMesh = refine(mesh)
#fineMesh = refine(fineMesh)
#fineMesh = refine(fineMesh)

fineSpace = FunctionSpace(fineMesh, "Lagrange", 1)
uF = interpolate(u, fineSpace)
eF = interpolate(e, fineSpace)

plot(uF)
plot(mesh)
interactive()

file = File('u.pvd')
file << uF
file = File('e.pvd')
file << eF
file = File('mesh.pvd')
file << mesh
