from dolfin import *

# Print log messages only from the root process in parallel
parameters["std_out_all_processes"] = False;

# Load mesh from file
N = 1
p_vel = 1
mesh = UnitSquareMesh(N,N)
#mesh = Mesh("lshape.xml.gz")

V = VectorFunctionSpace(mesh, "Quadrature", 1)
#V = FunctionSpace(mesh, "CG", p_vel)

# Define trial and test functions
u = TrialFunction(V)
v = TestFunction(V)

f = Expression(('sin(x[0]*pi)','cos(x[0]*pi)'))
#f = Expression('sin(x[0]*pi)')

A = assemble(inner(u,v)*dx) # mass matrix
b = assemble(inner(f,v)*dx)

print A.array()

# regular solve
u_h = Function(V)
solve(A, u_h.vector(), b)

# lumped solve
u_lumped = Function(V)
ones = Function(V)
ones.vector()[:] = 1.0
Mlumped = A*ones.vector()
u_lumped.vector().set_local(b.array()/Mlumped.array())

print u_h.vector().array() 
print u_lumped.vector().array()
print all(abs(u_h.vector().array() - u_lumped.vector().array()) < 1e-8)

