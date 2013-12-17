from dolfin import *
from numpy import array, argsort, amax
from math import ceil

useBulkChasing = True
useAdaptivity = False
numRefs = 4

# Create mesh and definene function space
pU = 1
pV = 4
N = 8
mesh = UnitSquareMesh(N,N)

# define problem params
eps = 1e-3
print "eps = ", eps
beta = Expression(('1.0','0.0'))
ue = Expression('(exp(r1*(x[0]-1))-exp(r2*(x[0]-1)))/(exp(-r1)-exp(-r2))*sin(pi*x[1])',eps = eps,pi = 2*acos(0.0), r1 = (-1+sqrt(1+4*pi*pi*eps*eps))/(-2*eps),r2 = (-1-sqrt(1+4*pi*pi*eps*eps))/(-2*eps))

zero = Expression('0.0')

def u0_boundary(x, on_boundary):
	return on_boundary
def outflow(x):
	return abs(x[0]-1) < 1E-14 or abs(x[1]-1) < 1E-14 or abs(x[1]) < 1E-14
def inflow(x):
	return outflow(x) == False;

class Inflow(Expression):
	def eval(self, values, x):
		values[0] = 0.0
		if abs(x[0]) < 1E-14:
			values[0] = 1.0
infl = Inflow()
outfl = 1-infl

hVec = []
errVec = []
nDofs = []
enrgy = []
for refIndex in xrange(numRefs):
	# define spaces
	U = FunctionSpace(mesh, "Lagrange", pU)
	V = FunctionSpace(mesh, "Lagrange", pV)
	E = U*V
	(u,e) = TrialFunctions(E)
	(du,v) = TestFunctions(E)
	n = FacetNormal(mesh)
	h = CellSize(mesh)

	bcs = []
	bcs.append(DirichletBC(E.sub(1), zero, outflow)) # error boundary condition
	bcs.append(DirichletBC(E.sub(0), ue, inflow)) # boundary conditions on u

	def ip(e,v):
		return inner(dot(beta,grad(e)),dot(beta,grad(v)))*dx + eps*inner(grad(e),grad(v))*dx #+ eps*inner(infl*e,v)*ds + eps*inner(outfl*dot(grad(e),n),dot(grad(v),n))*ds
	
	def b(u,v):
		return inner(-u,dot(beta,grad(v)))*dx + eps*inner(grad(u),grad(v))*dx - inner(eps*infl*dot(grad(u),n),v)*ds - inner(outfl*dot(beta,n)*u,v)*ds + inner(eps*outfl*u,dot(grad(v),n))*ds

	a = b(u,v) + b(du,e) + ip(e,v) #+ eps*inner(outfl*u,du)*ds

	f = Expression('0.0')
	x = V.cell().x
	L = inner(f,v)*dx - inner(infl*dot(beta,n)*ue,v)*ds
#	L = inner(f,v)*dx
	
	uSol = Function(E)
	solve(a==L, uSol, bcs)
	(u,e) = uSol.split()

	# evaluate element error indicator
	DG0 = FunctionSpace(mesh, "DG", 0) # element indicator function
	w = TestFunction(DG0)
	M = ip(w*e,e)
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
	
	#error = (u-ue)**2*dx
	#Err = sqrt(assemble(error))
	fineMesh = refine(mesh)
	fineSpace = FunctionSpace(fineMesh, "Lagrange", pU)
	uF = interpolate(u,fineSpace)
	Err = ue-uF;
	err = sqrt(assemble(Err**2*dx))
	h = (1.0/N)*.5**float(refIndex)
	
	print "on refinement ", refIndex
	energy = ip(e,e)
	totalE = sqrt(assemble(energy))
	print "h, ", h , ", E = ", err, ", e = ", totalE	
	hVec.append(h)
	errVec.append(err)
	nDofs.append(U.dofmap().global_dimension())
	enrgy.append(totalE)
	if (useAdaptivity):
		mesh = refine(mesh,cell_markers)
	else:
		mesh = refine(mesh)
	
# Convergence rates
from math import log as ln  # (log is a dolfin name too - and logg :-)
print 'h(', 1, ') = ', hVec[0], '; e(',1,') = ', errVec[0],'; energy(',1,') = ', enrgy[0], '; ndofs(',1,')=',nDofs[0]
for i in range(1, len(errVec)):
	r = ln(errVec[i]/errVec[i-1])/ln(hVec[i]/hVec[i-1])
	print 'h(', i+1, ') = ', hVec[i], '; e(',i+1,') = ', errVec[i],'; energy(',i+1,') = ', enrgy[i],'; r(',i,') = ', r, '; ndofs(',i+1,')=',nDofs[i]

for i in range(0, len(errVec)):
	print 'ratio(',i+1,') = ', errVec[i]/(enrgy[i])


# Plot solution
fineMesh = refine(mesh)
#fineMesh = refine(fineMesh)
#fineMesh = refine(fineMesh)
fineSpace = FunctionSpace(fineMesh, "Lagrange", 1)
uF = interpolate(u, fineSpace)
eF = interpolate(e, fineSpace)
err = ue-uF

plot(uF)
#plot(eF)
plot(err)
interactive()


#file = File('u.pvd')
#file << uF
#file = File('e.pvd')
#file << eF
#file = File('mesh.pvd')
#file << mesh
