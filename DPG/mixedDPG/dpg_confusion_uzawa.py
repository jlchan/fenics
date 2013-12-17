from dolfin import *

# Create mesh and define function space
pU = 1
pV = 3
N = 16
mesh = UnitSquareMesh(N,N)

# define problem params
beta = Expression(('1.0','0.0'))
eps = 1e-3
uex = Expression('sin(x[1]*pi*2)')
#uex = Expression('(exp(r1*(x[0]-1))-exp(r2*(x[0]-1)))/(exp(-r1)-exp(-r2))*sin(pi*x[1])',eps = eps,pi = 2*acos(0.0), r1 = (-1+sqrt(1+4*pi*pi*eps*eps))/(-2*eps),r2 = (-1-sqrt(1+4*pi*pi*eps*eps))/(-2*eps))

zero = Constant('0.0')

def u0_boundary(x, on_boundary):
	return on_boundary
def inflow(x):
	return abs(x[0]) < 1E-14 
def outflow(x):
	return abs(x[0]-1) < 1E-14 or abs(x[1]) < 1E-14 or abs(x[1]-1) < 1E-14
class Inflow(Expression):
	def eval(self, values, x):
		values[0] = 0.0
		if inflow(x):
			values[0] = 1.0
infl = Inflow()
outfl = Expression('1.0')-infl

# define spaces
U = FunctionSpace(mesh, "Lagrange", pU)
V = FunctionSpace(mesh, "Lagrange", pV)
E = U*V
(u,e) = TrialFunctions(E)
(du,v) = TestFunctions(E)
n = FacetNormal(mesh)
h = CellSize(mesh)

def ip(e,v):
	return inner(dot(beta,grad(e)),dot(beta,grad(v)))*dx + inner(e,v)*dx + eps*inner(grad(e),grad(v))*dx #+ inner(infl*e,v)*ds

def b(u,v):
	return inner(outfl*dot(beta,n)*u,v)*ds + inner(u,-dot(beta,grad(v)))*dx + eps*inner(grad(u),grad(v))*dx - inner(eps*infl*dot(grad(u),n),v)*ds

def l(v):
	f = Expression('0.0')
	return inner(f,v)*dx - inner(infl*dot(beta,n)*uex,v)*ds

def m(u,du):
	return inner(u,du)*dx + eps*inner(grad(u),grad(du))*dx + eps*inner(infl*dot(grad(u),n),dot(grad(duh),n))*ds

a = b(u,v) + b(du,e) + ip(e,v)
uSol = Function(E)
bcs = [DirichletBC(E.sub(0), uex, inflow), DirichletBC(E.sub(1), zero, outflow)]
#bcs = [DirichletBC(E.sub(0), uex, u0_boundary), DirichletBC(E.sub(1), zero, u0_boundary)]
solve(a==l(v), uSol, bcs)
(u,e) = uSol.split()

fineMesh = refine(mesh)
#fineMesh = refine(fineMesh)
#fineMesh = refine(fineMesh)

fineSpace = FunctionSpace(fineMesh, "Lagrange", pU)
uF = interpolate(u, fineSpace)
eF = interpolate(e, fineSpace)
err = uex-uF
l2err = sqrt(abs(assemble(err**2*dx)))
print "err = ", l2err
energyErr= sqrt(abs(assemble(ip(e,e))))
print "energy err = ", energyErr
print "l2/energy ratio = ", l2err/energyErr

# ===================== Trying Uzawa =========================

#(e,z)_V = <f-Au,z>
#(uh,v) = (uh,v) + (e,v)
r = TrialFunction(V)
z = TestFunction(V)
uh = TrialFunction(U)
duh = TestFunction(U)
rk = Function(V)
uhk = Function(U)
uhkPrev = Function(U)
uhk = interpolate(Expression('1.0'),U)
uhkPrev = interpolate(uhk,U)
bcU = DirichletBC(U, uex, inflow) # error boundary condition
bcV = DirichletBC(V, zero, outflow) # error boundary condition
#bcU = DirichletBC(U, uex, u0_boundary) # error boundary condition
#bcV = DirichletBC(V, zero, u0_boundary) # error boundary condition

relres = []
for i in range(0,25):
	R_V = ip(r,z)
	residual = l(z) - b(uhk,z)
	solve(R_V==residual,rk,bcV)
	M = inner(uh,duh)*dx	
	uzawa = inner(uhk,duh)*dx + b(duh,rk) #inner(Adjoint(rk),duh)*dx
	uhkPrev = interpolate(uhk,U)
	solve(M==uzawa,uhk,bcU)	
	duhk = sqrt(abs(assemble((uhk-uhkPrev)**2*dx)))
	print "on iter ", i, ", duh_k = ", duhk
	relres.append(duhk)

for i in range(1, len(relres)):
	print relres[i]

#exit(0)
# ===================== Trying Uzawa =========================

plot(interpolate(uhk,fineSpace))
plot(uF)
plot(err)
interactive()

#file = File('u.pvd')
#file << uF
#file = File('e.pvd')
#file << eF
#file = File('mesh.pvd')
#file << mesh
