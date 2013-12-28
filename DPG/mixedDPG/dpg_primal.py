from dolfin import *
from numpy import array, argsort, amax
from math import ceil

import helper_functions as help
import helper_confusion as chelp

useBulkChasing = False
useAdaptivity = False

# Create mesh and definene function space
eps = float(help.parseArg('--eps',1e-2))
pU = int(help.parseArg('--p',1))
N = int(help.parseArg('--N',4))
numRefs = int(help.parseArg('--numRefs',1))
numQRefs = int(help.parseArg('--numQRefs',4))

useStrongBC = help.parseArg('--useStrongBC','False')=='True' #eval using strings
plotFlag = help.parseArg('--plot','True')=='True' #eval using strings

dp = int(help.parseArg('--dp',1))

pV = pU+dp
mesh = UnitSquareMesh(N,N)

zero = Expression('0.0')

# define problem params
print "eps = ", eps
beta = Constant((.5,1.0))
ue = zero

def outflow(x):
	return abs(x[0]-1) < 1E-14 or abs(x[1]-1) < 1E-14 
def inflow(x):
	return abs(x[0]) < 1E-14 or abs(x[1]) < 1E-14 #or abs(x[1]-1) < 1E-14 or abs(x[1]) < 1E-14
class Inflow(Expression):
	def eval(self, values, x):
		values[0] = 0.0
		if inflow(x):
			values[0] = 1.0
infl = Inflow()
outfl = 1.0-infl

enrgy = []
for refIndex in xrange(numRefs):
	# define spaces
	U = FunctionSpace(mesh, "CG", pU)
	V = FunctionSpace(mesh, "DG", pV)
	if (pU>1):	
		F = FunctionSpace(mesh,"BDM",pU-1,restriction="facet") # interfacial flux		
	else:
		F = FunctionSpace(mesh,"RT",pU,restriction="facet") # interfacial flux

	E = MixedFunctionSpace( [U,V,F] )
	(u,e,f) = TrialFunctions(E)
	(du,v,df) = TestFunctions(E)
	n = FacetNormal(mesh)
	h = CellSize(mesh)

	bcs = []
	bcs.append(DirichletBC(E.sub(2), Constant((0.0,0.0)), DomainBoundary())) # flux BCs - remove their contribution completely!
	bcs.append(DirichletBC(E.sub(1), zero, outflow,"pointwise")) # error boundary condition
	bcs.append(DirichletBC(E.sub(0), ue, inflow)) # boundary conditions on u
	if useStrongBC:
		bcs.append(DirichletBC(E.sub(0), ue, outflow)) # boundary conditions on u

	def ip(e,v):
		return inner(dot(beta,grad(e)),dot(beta,grad(v)))*dx + eps*inner(grad(e),grad(v))*dx + inner(e,v)*dx	
	def b(u,v,f):
		fieldForm = inner(dot(beta,n)*u,v)*ds + inner(-u,dot(beta,grad(v)))*dx + eps*inner(grad(u),grad(v))*dx 
		fieldForm += inner(dot(f('+'),n('+')),v('+')-v('-')) * dS 
		if useStrongBC: # strong outflow
			return fieldForm - inner(eps*dot(grad(u),n),v)*ds 
		else: # nitsche type of weak BC
			return fieldForm - inner(eps*dot(grad(u),n),v)*ds - inner(eps*outfl*u,dot(grad(v),n))*ds
	def pen(f,df):
		bn = abs(dot(beta('+'),n('+')))
		alpha = 1.0/(bn+eps + 1e-7)
		return .001*inner(alpha*dot(f('+'),n('+')),dot(df('+'),n('+'))) * dS 
		#return inner(alpha*f('+'),df('+')) * dS 

	a = b(u,v,f) + b(du,e,df) + ip(e,v) + pen(f,df)

	x = V.cell().x
	forcing = conditional(le( x[1]-2.0*x[0],  0.0), 1.0, 0.0) # discontinuous forcing data
	L = inner(forcing,v)*dx 
	
	uSol = Function(E)
	solve(a==L, uSol, bcs)
	(uh,eh,fh) = uSol.split()

	# evaluate element error indicator
	DG0 = FunctionSpace(mesh, "DG", 0) # element indicator function
	w = TestFunction(DG0)
	cell_energy = assemble(ip(w*eh,eh))	

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
	
	print "on refinement ", refIndex
	energy = ip(eh,eh)
	totalE = sqrt(abs(assemble(energy)))

	enrgy.append(totalE)
	if (useAdaptivity):
		mesh = refine(mesh,cell_markers)
	else:
		mesh = refine(mesh)
	
print "energy_err = ", enrgy

if plotFlag:
	#print "solving a final time"
	#uSol = Function(E)
	#solve(a==L, uSol, bcs)
	#(uh,eh,fh) = uSol.split()

	# Plot solution
	fineMesh = mesh
	for ref in xrange(pU-1):
		fineMesh = refine(fineMesh)
	fineSpace = FunctionSpace(fineMesh, "Lagrange", 1)
	eF = interpolate(eh, fineSpace)
	uF = interpolate(uh, fineSpace)
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

