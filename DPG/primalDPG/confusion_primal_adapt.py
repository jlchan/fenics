from dolfin import *

numRefs = 4

# Create mesh and define function space
degree = 2
delta_p = 2
N = 8
mesh = UnitSquareMesh(N,N)

# define problem params
eps = 1.0e-4
print "eps = ", eps
beta = Constant((1.0,1.0))#Constant((.5,1.0))
zero = Constant(0.0)

u0 = zero #Expression('sin(x[0]*pi*4)')
# restricts to inflow only
class InflowInd(Expression):
	def eval(self, values, x):
		values[0] = 0.0
		if abs(x[0]) < 1E-14 or abs(x[1]) < 1E-14:
			values[0] = 1.0
inflowIndicator = InflowInd()
outflowIndicator = Expression('1.0')-inflowIndicator

class Inflow(SubDomain):	def inside(self,x, on_boundary):		tol = 1E-14		return on_boundary and (abs(x[0]) < tol or abs(x[1]) < tol)
class Outflow(SubDomain):	def inside(self,x, on_boundary):		tol = 1E-14		return on_boundary and (abs(x[0]-1.0) < tol or abs(x[1]-1.0) < tol)

inflow = Inflow()
outflow = Outflow()
boundary_parts = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
inflow.mark(boundary_parts,0)
outflow.mark(boundary_parts,1)

for refIndex in xrange(numRefs):
	# define spaces
	U = FunctionSpace(mesh, "CG", degree+1)					# trial functions
	V = FunctionSpace(mesh, "DG", degree+delta_p)						# discontinuous test functions
	F = FunctionSpace(mesh, "BDM", degree, restriction="facet") 	# fluxes
	E = MixedFunctionSpace([U,V,F])
	(u,e,f) = TrialFunctions(E)
	(du,v,df) = TestFunctions(E)
	n = FacetNormal(mesh)
	h = CellSize(mesh)
	
	def ip(e,v):
		w = Expression('x[0]*x[1]') + Constant(eps)
		return w*inner(dot(beta,grad(e)),dot(beta,grad(v))) + eps*inner(grad(e),grad(v)) + inner(e,v)
	
	def b(u,v,f):
		b = inner(-u,dot(beta,grad(v)))*dx + eps*inner(grad(u),grad(v))*dx 
		b = b - inner(dot(f('+'),n('+')),jump(v)) * dS - inner(dot(f,n),v)*ds
		return b
	
	a = b(u,v,f) + b(du,e,df) + ip(e,v)*dx

	load = Constant(1.0) #zero
	L = inner(load,v)*dx
#	L = L - inner(inflowIndicator*dot(beta,n)*u0,v)*ds
	
	bcs = []
	bcs.append(DirichletBC(E.sub(0), u0, DomainBoundary())) # boundary conditions on u
#	bcs.append(DirichletBC(E.sub(2), Constant((0.0,0.0)), boundary_parts, 0)) # boundary conditions on u
#	bcs.append(DirichletBC(E.sub(0), zero, outflow)) # bcs on error
#	bcs.append(DirichletBC(E.sub(1), zero, outflow)) # bcs on error
	

#	solve(a==L, uSol, bcs)
	A = assemble(a, exterior_facet_domains=boundary_parts)	b = assemble(L, exterior_facet_domains=boundary_parts)
	for bc in bcs: bc.apply(A,b)
	uSol = Function(E)
	Uv = uSol.vector()
#	solve(a==L, uSol)
	solve(A,Uv,b)
	(u,e,f) = uSol.split()
	
	# evaluate element error indicator
	DG0 = FunctionSpace(mesh, "DG", 0) # element indicator function
	I_K = TestFunction(DG0)
	M = ip(e,e)*I_K*dx 
	cell_energy = assemble(M)	
	E = sqrt(assemble(ip(e,e)*dx))
	num_trial_dofs = U.dofmap().global_dimension() + F.dofmap().global_dimension()

	# define adaptive strategy
	factor = .25
	cell_markers = MeshFunction("bool", mesh, mesh.topology().dim())
	cell_markers.set_all(False)
	cutoff = sorted(cell_energy, reverse=True)[int(len(cell_energy)*factor)]
	for c in cells(mesh):
		cell_markers[c] = cell_energy[c.index()] > cutoff
	
	mesh = refine(mesh,cell_markers)
#	mesh = refine(mesh)
	print "on refinement ", refIndex, " energy = ", E, " with num trial dofs = ", num_trial_dofs
		
fineMesh = refine(refine(mesh))
fineSpace = FunctionSpace(fineMesh, "Lagrange", 1)
uF = interpolate(u, fineSpace)
eF = interpolate(e, fineSpace)

plot(uF)
plot(eF)
plot(mesh)
interactive()

#file = File('u.pvd')
#file << uF
#file = File('e.pvd')
#file << eF
#file = File('mesh.pvd')
#file << mesh
