from dolfin import *
from numpy import array, argsort, amax
from math import ceil

import helper_functions as help
import helper_confusion as chelp

useBulkChasing = True
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

# define problem params
print "eps = ", eps
beta = chelp.beta()
ue = chelp.erikkson_solution(eps)
grad_ue = chelp.erikkson_solution_gradient(eps)

zero = Expression('0.0')

infl = chelp.Inflow()
outfl = 1-infl

hVec = []
errVec = []
nitscheErrVec = []
nDofs = []
enrgy = []
for refIndex in xrange(numRefs):
	# define spaces
	U = FunctionSpace(mesh, "Lagrange", pU)
	V = FunctionSpace(mesh, "DG", pV)
	F = FunctionSpace(mesh,"RT",pU,restriction="facet") # interfacial flux q                                                 
	E = MixedFunctionSpace( [U,V,F] )
	(u,e,f) = TrialFunctions(E)
	(du,v,df) = TestFunctions(E)
	n = FacetNormal(mesh)
	h = CellSize(mesh)

	bcs = []
	bcs.append(DirichletBC(E.sub(1), zero, chelp.outflow,"pointwise")) # error boundary condition
	bcs.append(DirichletBC(E.sub(0), ue, chelp.inflow)) # boundary conditions on u
	if useStrongBC:
		bcs.append(DirichletBC(E.sub(0), ue, chelp.outflow)) # boundary conditions on u

	def ip(e,v):
		return inner(dot(beta,grad(e)),dot(beta,grad(v)))*dx + eps*inner(grad(e),grad(v))*dx + inner(e,v)*dx	
	def b(u,v,f):
		fieldForm = inner(dot(beta,n)*u,v)*ds + inner(-u,dot(beta,grad(v)))*dx + eps*inner(grad(u),grad(v))*dx 
		fluxForm = -inner(dot(f('+'),n('+')),v('+')-v('-')) * dS - inner(dot(f,n),v) * ds
		if useStrongBC: # strong outflow
			return fieldForm + fluxForm - inner(eps*dot(grad(u),n),v)*ds 
		else: # nitsche type of weak BC
			return fieldForm + fluxForm - inner(eps*dot(grad(u),n),v)*ds - inner(eps*outfl*u,dot(grad(v),n))*ds

	a = b(u,v,f) + b(du,e,df) + ip(e,v) #+ inner(dot(f,n),dot(df,n))*ds

	F = Expression('0.0')
	x = V.cell().x
	L = inner(F,v)*dx 
	
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
	
	fineMesh = chelp.quadrature_refine(mesh,N,numQRefs)
	fineSpace = FunctionSpace(fineMesh, "Lagrange", pU+2)
	uF = interpolate(uh,fineSpace)
	L2err = (ue-uF)**2*dx
	l_err = sqrt(abs(assemble(L2err)))
	hh = (1.0/N)*.5**float(refIndex)
	
	print "on refinement ", refIndex
	energy = ip(eh,eh)
	totalE = sqrt(abs(assemble(energy)))

	H1err = (grad_ue-grad(uF))**2*dx
	edge_error = infl*1/h*(ue-uF)**2*ds
	n_err = sqrt(abs(assemble(eps*edge_error + eps*H1err + L2err)))
	nitscheErrVec.append(n_err)

	print "h, ", hh , ", L2 error = ", l_err, ", e = ", totalE	
	hVec.append(hh)
	errVec.append(l_err)
	nDofs.append(U.dofmap().global_dimension())
	enrgy.append(totalE)
	if (useAdaptivity):
		mesh = refine(mesh,cell_markers)
	else:
		mesh = refine(mesh)
	
# Convergence rates
from math import log as ln  # (log is a dolfin name too - and logg :-)
#print 'h(', 1, ') = ', hVec[0], '; e(',1,') = ', errVec[0],'; energy(',1,') = ', enrgy[0], '; ndofs(',1,')=',nDofs[0]
rVec = []
for i in range(1, len(errVec)):
	r = ln(errVec[i]/errVec[i-1])/ln(hVec[i]/hVec[i-1])
	rVec.append(r)
	#print 'h(', i+1, ') = ', hVec[i], '; e(',i+1,') = ', errVec[i],'; energy(',i+1,') = ', enrgy[i],'; r(',i,') = ', r, '; ndofs(',i+1,')=',nDofs[i]

ratVec = []
nratVec = []
for i in range(0, len(errVec)):
	ratVec.append(errVec[i]/enrgy[i])
	nratVec.append(nitscheErrVec[i]/enrgy[i])
	#print 'ratio(',i+1,') = ', errVec[i]/(enrgy[i])

print "h = ", hVec
print "energy_err = ", enrgy
print "L2_err = ", errVec
print "rates = ", rVec
print "L2_ratios = ", ratVec
print "nitsche_ratios = ", nratVec

if plotFlag:
	# Plot solution
	fineMesh = mesh
	for ref in xrange(pU-1):
		fineMesh = refine(fineMesh)
		fineSpace = FunctionSpace(fineMesh, "Lagrange", 1)
	uF = interpolate(uh, fineSpace)
	eF = interpolate(eh, fineSpace)
	err = ue-uF
	plot(uF)
	plot(eF)
	plot(mesh)
	#plot(err)
	interactive()


#file = File('u.pvd')
#file << uF
#file = File('e.pvd')
#file << eF
#file = File('mesh.pvd')
#file << mesh
