
from dolfin import *

# Print log messages only from the root process in parallel
parameters["std_out_all_processes"] = False;

# Load mesh from file
N = 16
mesh = UnitSquareMesh(N,N)

# Define function spaces (P2-P1)
V = VectorFunctionSpace(mesh, "CG", 2)
Q = FunctionSpace(mesh, "CG", 1)

# Define trial and test functions
u = TrialFunction(V)
p = TrialFunction(Q)
v = TestFunction(V)
q = TestFunction(Q)

# Set parameter values
dt = 0.01
T = 10
nu = 1.0/40.0

# Define time-dependent pressure boundary condition
#p_in = Expression("sin(3.141592*t)", t=0.0)
Re2 = 2/nu
lam = Re2 - sqrt(Re2**2 + 4*pi**2)
ue = Expression(("1-exp(lam*x[0])*cos(2*pi*x[1])","lam/(2*pi)*exp(lam*x[0])*sin(2*pi*x[1])"),pi = pi,lam = lam)
pe_no_avg = Expression(".5*exp(2*lam*x[0])", lam=lam)
p_avg = Expression("(exp(2.0*lam) - 1.0)/(4.0*lam)",lam=lam)
pe = pe_no_avg #- p_avg
ubc = DirichletBC(V, ue,"on_boundary")
pbc = DirichletBC(Q, pe,"x[0] > 1.0 - DOLFIN_EPS | x[0] < DOLFIN_EPS")
bcu = [ubc]
bcp = [pbc]

# Create functions to store intermediate values
u0 = Function(V)
u0 = interpolate(ue,V)
u1 = Function(V)
p1 = Function(Q)

# Define coefficients
k = Constant(dt)
f = Constant((0, 0))

# Tentative velocity step
F1 = (1/k)*inner(u - u0, v)*dx + inner(grad(u0)*u0, v)*dx + \
     nu*inner(grad(u), grad(v))*dx - inner(f, v)*dx
a1 = lhs(F1)
L1 = rhs(F1)

# Pressure update
a2 = inner(grad(p), grad(q))*dx
L2 = -(1/k)*div(u1)*q*dx

# Velocity update
a3 = inner(u, v)*dx
L3 = inner(u1, v)*dx - k*inner(grad(p1), v)*dx

# Assemble matrices
A1 = assemble(a1)
A2 = assemble(a2)
A3 = assemble(a3)

# Use amg preconditioner if available
prec = "amg" if has_krylov_solver_preconditioner("amg") else "default"

# Create files for storing solution
ufile = File("results/velocity.pvd")
pfile = File("results/pressure.pvd")

plot(interpolate(pe,Q), title="Exact Pressure",rescale=True) 
plot(interpolate(ue,V), title="Exact Velocity",rescale=True) 

# Time-stepping
t = dt
while t < T + DOLFIN_EPS:

    # Update pressure boundary condition
    #p_in.t = t

    # Compute tentative velocity step
    begin("Computing tentative velocity")
    b1 = assemble(L1)
    [bc.apply(A1, b1) for bc in bcu]
    solve(A1, u1.vector(), b1, "gmres", "default")
    end()

    # Pressure correction
    begin("Computing pressure correction")
    b2 = assemble(L2)
    [bc.apply(A2, b2) for bc in bcp]
    solve(A2, p1.vector(), b2, "gmres", prec)
    #p1.assign(interpolate(pe,Q))
    end()

    # Velocity correction
    begin("Computing velocity correction")
    b3 = assemble(L3)
    [bc.apply(A3, b3) for bc in bcu]
    solve(A3, u1.vector(), b3, "gmres", "default")
    end()

    # Plot solution
    plot(p1, title="Pressure", rescale=True)
    plot(u1, title="Velocity", rescale=True)

    # Save to file
    ufile << u1
    pfile << p1

    # Move to next time step
    u0.assign(u1)
    t += dt
    print "t =", t

# Hold plot
interactive()
