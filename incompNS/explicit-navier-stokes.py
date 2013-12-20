from dolfin import *

# Print log messages only from the root process in parallel
parameters["std_out_all_processes"] = False;

# Load mesh from file
N = 32
p_vel = 2
mesh = UnitSquareMesh(N,N)
hh = 1.0/(2.0*N) #circrumscribed circle defn

# Define function spaces (P2-P1)
V = VectorFunctionSpace(mesh, "CG", p_vel)
Q = FunctionSpace(mesh, "CG", 1)

# Define trial and test functions
u = TrialFunction(V)
p = TrialFunction(Q)
v = TestFunction(V)
q = TestFunction(Q)

# Set parameter values
dt = .01#.1*hh/(p_vel**2)
print "dt = ", dt
T = 100
nu = 1.0/40.0

# Define time-dependent pressure boundary condition
#p_in = Expression("sin(3.141592*t)", t=0.0)
#p_out = Expression("0.0")
# Define boundary conditions
#noslip  = DirichletBC(V, (0, 0),"on_boundary && (x[0] < DOLFIN_EPS | x[0] > 1.0 - DOLFIN_EPS )")
#inflow  = DirichletBC(Q, p_in, "x[1] > 1.0 - DOLFIN_EPS")
#outflow = DirichletBC(Q, p_out, "x[1] < DOLFIN_EPS")
#bcu = [noslip]
#bcp = [inflow, outflow]
Re2 = 2/nu
lam = Re2 - sqrt(Re2**2 + 4*pi**2)
ue = Expression(("1-exp(lam*x[0])*cos(2*pi*x[1])","lam/(2*pi)*exp(lam*x[0])*sin(2*pi*x[1])"),pi = pi,lam = lam)
pe_no_avg = Expression(".5*exp(2*lam*x[0])", lam=lam)
p_avg = Expression("(exp(2.0*lam) - 1.0)/(4.0*lam)",lam=lam)
pe = pe_no_avg #- p_avg
p_in = Expression("1.0")
ubc=DirichletBC(V, ue,"on_boundary")
pbc=DirichletBC(Q, pe,"x[0] > 1.0 - DOLFIN_EPS | x[0] < DOLFIN_EPS")
bcu = [ubc]
bcp = [pbc]

# Create functions to store intermediate values
u0 = Function(V)
u0 = interpolate(ue,V)
u1 = Function(V)
u2 = Function(V)
p1 = Function(Q)

# Define coefficients
k = Constant(dt)
f = Constant((0, 0))

# Tentative velocity step
F1 = (1/k)*inner(u - u0, v)*dx + inner(grad(u0)*u0, v)*dx - inner(f, v)*dx
a1 = lhs(F1)
L1 = rhs(F1)

# Pressure update
a2 = inner(grad(p), grad(q))*dx
L2 = -(1/k)*div(u1)*q*dx

# Velocity update
a3 = inner(u, v)*dx
L3 = inner(u1, v)*dx - k*inner(grad(p1), v)*dx

# Viscous update
a4 = inner(grad(u), grad(v))*dx + (1/(k*nu))*inner(u,v)*dx
L4 = (1/(k*nu))*inner(u1, v)*dx

# Assemble matrices
A1 = assemble(a1)
A2 = assemble(a2)
A3 = assemble(a3)
A4 = assemble(a4)

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
    p_in.t = t

    # Compute tentative velocity step
    begin("Convecting explicitly velocity")
    b1 = assemble(L1)
    [bc.apply(A1, b1) for bc in bcu]
    solve(A1, u1.vector(), b1, "gmres", "default")
    #lumped = assemble(action(a1), coefficients=[Constant(1)])
    #Uex = u1.vector()
    #for i, bi in enumerate(b1):
    #    Uex[i] = bi/lumped[i]
    #u1.vector().set_local(Uex[V.dofmap().vertex_to_dof_map(mesh)])
    end()

    # Pressure correction
    begin("Computing pressure correction")
    b2 = assemble(L2)
    [bc.apply(A2, b2) for bc in bcp]
    solve(A2, p1.vector(), b2, "gmres", prec)
    end()

    # Velocity correction
    begin("Computing velocity correction")
    b3 = assemble(L3)
    [bc.apply(A3, b3) for bc in bcu]
    solve(A3, u1.vector(), b3, "gmres", "default")
    end()

    begin("Computing viscous correction")
    b4 = assemble(L4)
    [bc.apply(A4, b4) for bc in bcu]
    solve(A4, u2.vector(), b4, "gmres", "default")
    end()

    # Plot solution
    plot(p1, title="Pressure", rescale=True)
    plot(u2, title="Velocity", rescale=True)

    # Save to file
    ufile << u2
    pfile << p1

    # Move to next time step
    u0.assign(u2)
    t += dt
    print "t =", t

# Hold plot
interactive()
