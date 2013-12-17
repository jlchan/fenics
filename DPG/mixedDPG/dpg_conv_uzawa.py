from dolfin import *

# Create mesh and define function space
pU = 1
pV = 2
N = 8
mesh = UnitSquareMesh(N,N)

# define problem params
beta = Expression(('1.0','0.0'))

uex = Expression('sin(x[1]*pi*2)')
zero = Constant('0.0')

def u0_boundary(x, on_boundary):
	return on_boundary
def inflow(x):
	return abs(x[0]) < 1E-14 
def outflow(x):
	return abs(x[1]) < 1E-14
class Inflow(Expression):
	def eval(self, values, x):
		values[0] = 0.0
		if (abs(x[0]) < 1E-14 or abs(x[0]-1) < 1E-14 or abs(x[1]-1) < 1E-14): #or abs(x[1]) < 1E-14 or abs(x[1]-1.0) < 1E-14:
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
	return inner(dot(beta,grad(e)),dot(beta,grad(v)))*dx + inner(e,v)*dx #+ 1e-1*inner(grad(e),grad(v))*dx

def b(u,v):
	return inner(outfl*dot(beta,n)*u,v)*ds + inner(u,-dot(beta,grad(v)))*dx

def l(v):
	f = Expression('0.0')
	return inner(f,v)*dx - inner(infl*dot(beta,n)*uex,v)*ds

def Adjoint(r):
	return -dot(beta,grad(r))

a = b(u,v) + b(du,e) + ip(e,v)
uSol = Function(E)
solve(a==l(v), uSol)
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
for i in range(0,10):
	R_V = ip(r,z)
	residual = l(z) - b(uhk,z)
	solve(R_V==residual,rk)
	M = inner(uh,duh)*dx	
	uzawa = inner(uhk,duh)*dx + inner(Adjoint(rk),duh)*dx
	uhkPrev = interpolate(uhk,U)
	solve(M==uzawa,uhk)	
	duhk = sqrt(abs(assemble((uhk-uhkPrev)**2*dx)))
	print "on iter ", i, ", duh_k = ", duhk

# ===================== Trying Uzawa =========================
plot(interpolate(uhk,fineSpace))
plot(uF)
plot(err)
interactive()

file = File('u.pvd')
file << uF
file = File('e.pvd')
file << eF
file = File('mesh.pvd')
file << mesh
