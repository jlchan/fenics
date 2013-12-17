from dolfin import *

# Create mesh and define function space
pU = 2
pV = 3
N = 8
mesh = UnitSquareMesh(N,N)

# define problem params
beta = Expression(('1.0','0.0'))
eps = 1e-0
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
x = V.cell().x

def ip(e,v):
	return inner(dot(beta,grad(e)),dot(beta,grad(v)))*dx + eps*inner(e,v)*dx + eps*inner(grad(e),grad(v))*dx + inner(infl*e,v)*ds

def b(u,v):
	return inner(outfl*dot(beta,n)*u,v)*ds + inner(u,-dot(beta,grad(v)))*dx + eps*inner(grad(u),grad(v))*dx - inner(eps*infl*dot(grad(u),n),v)*ds

def l(v):
	f = Expression('0.0')
	return inner(f,v)*dx - inner(infl*dot(beta,n)*uex,v)*ds

def Q(du):
	f = conditional(le( (x[1]-.5)**2 + (x[0]-.5)**2,  .01), 1.0, 0.0) # discontinuous forcing data in circle
	return inner(f,du)*dx


def m(u,du):
	return inner(u,du)*dx + inner(dot(grad(u),beta),dot(grad(du),beta))*dx + eps*inner(grad(u),grad(du))*dx + eps*inner(infl*dot(grad(u),n),dot(grad(duh),n))*ds

a = b(u,v) + b(du,e) + ip(e,v)
uSol = Function(E)
bcs = [DirichletBC(E.sub(0), zero, inflow), DirichletBC(E.sub(1), zero, outflow)]
#bcs = [DirichletBC(E.sub(0), zero, u0_boundary), DirichletBC(E.sub(1), zero, u0_boundary)]
solve(a==Q(-du), uSol, bcs)
(u,e) = uSol.split()

fineMesh = refine(mesh)
fineMesh = refine(fineMesh)
#fineMesh = refine(fineMesh)

fineSpace = FunctionSpace(fineMesh, "Lagrange", pU)
uF = interpolate(u, fineSpace)
eF = interpolate(e, fineSpace)

# ===================== Trying Uzawa =========================

#(e,z)_V = <f-Au,z>
#(uh,v) = (uh,v) + (e,v)
r = TrialFunction(V)
z = TestFunction(V)
uh = TrialFunction(U)
duh = TestFunction(U)
rk = Function(V)
uhk = Function(U)
#bcU = DirichletBC(U, zero, inflow) # error boundary condition
#bcV = DirichletBC(V, zero, outflow) # error boundary condition
bcU = DirichletBC(U, zero, u0_boundary) # error boundary condition
bcV = DirichletBC(V, zero, u0_boundary) # error boundary condition

M = m(uh,duh)
b = Q(duh) #inner(uhk,duh)*dx + b(duh,rk) #inner(Adjoint(rk),duh)*dx
solve(M==b,uhk,bcU)	
#solve(M==b,uhk)	

#exit(0)
# ===================== Trying Uzawa =========================

plot(uF)
plot(interpolate(uhk,fineSpace))

interactive()

#file = File('u.pvd')
#file << uF
#file = File('e.pvd')
#file << eF
#file = File('mesh.pvd')
#file << mesh
