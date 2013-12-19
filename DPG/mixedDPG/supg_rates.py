from dolfin import *
import helper_functions as help

# define problem params
eps = 1e-4
print "eps = ", eps
beta = Constant(('1.0','0.0'))
ue = Expression('(exp(r1*(x[0]-1))-exp(r2*(x[0]-1)))/(exp(-r1)-exp(-r2))*sin(pi*x[1])',eps = eps,pi = 2*acos(0.0), r1 = (-1+sqrt(1+4*pi*pi*eps*eps))/(-2*eps),r2 = (-1-sqrt(1+4*pi*pi*eps*eps))/(-2*eps))

zero = Expression('0.0')

def u0_boundary(x, on_boundary):
	return on_boundary

N = 4
mesh = UnitSquareMesh(N,N)
numRefs = 5

# Load mesh and subdomains
p = 2
for refIndex in xrange(numRefs):
    h = CellSize(mesh)
    
    # Create FunctionSpaces
    V = FunctionSpace(mesh, "CG", p)
    
    # Test and trial functions
    u = TrialFunction(V) 
    v = TestFunction(V)
    
    f = Constant(0.0)
    
    # Galerkin variational problem
    F = inner(u,-dot(beta,grad(v)))*dx + eps*inner(grad(u),grad(v))*dx - inner(f,v)*dx
    
    # Add SUPG stabilisation terms
    bnorm = sqrt(dot(beta,beta))
    residual = div(beta*u) - eps*div(grad(u)) - f
    F += h/((p+1)*bnorm)*dot(beta, grad(v))*residual*dx
    
    # Create bilinear and linear forms
    a = lhs(F)
    L = rhs(F)
    
    bc = DirichletBC(V, ue, u0_boundary)

    uh = Function(V)
    solve(a==L,uh,bc)

    # compute error on refined mesh for accuracy
    tempMesh = mesh
    for i in xrange(p):
        tempMesh = refine(tempMesh)

    Vp = FunctionSpace(tempMesh,"CG",p+1)
    uh = interpolate(uh,Vp)
    error = (ue-uh)**2*dx
    e = sqrt(assemble(error))
    hh = (1.0/N)*.5**float(refIndex)

    print "for h = ", hh, ", L2 error = ", e
    
    # refine mesh for next iteration
    mesh = refine(mesh)

plot(uh)
interactive()
