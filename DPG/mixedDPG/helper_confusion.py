from dolfin import * 

# =======================   FEM SPECIFIC ====================================

# refines a mesh at the inflow to improve quadrature
def quadrature_refine(mesh,N,nqrefs):
    tempMesh = mesh
    for nqr in xrange(nqrefs):
        in_cells = MeshFunction("bool", tempMesh, tempMesh.topology().dim())
        in_cells.set_all(False)
        for c in cells(tempMesh):
            in_cells[c] = any([near(vertex.x(0),1) for vertex in vertices(c)]) 
        tempMesh = refine(tempMesh,in_cells)
    return tempMesh

def erikkson_solution(eps):
    r1=(-1+sqrt(1+4*pi*pi*eps*eps))/(-2*eps)
    r2=(-1-sqrt(1+4*pi*pi*eps*eps))/(-2*eps)
    return Expression('(exp(r1*(x[0]-1))-exp(r2*(x[0]-1)))*C*sin(pi*x[1])',eps = eps,pi = 2*acos(0.0), r1=r1, r2=r2, C = 1/(exp(-r1)-exp(-r2)))

def erikkson_solution_gradient(eps):
    r1=(-1+sqrt(1+4*pi*pi*eps*eps))/(-2*eps)
    r2=(-1-sqrt(1+4*pi*pi*eps*eps))/(-2*eps)
    C = 1/(exp(-r1)-exp(-r2))
    ux = '(r1*exp(r1*(x[0]-1))-r2*exp(r2*(x[0]-1)))*C*sin(pi*x[1])'
    uy = '(exp(r1*(x[0]-1))-exp(r2*(x[0]-1)))*C*pi*cos(pi*x[1])'
    return Expression((ux,uy),eps = eps,pi = 2*acos(0.0), r1=r1, r2=r2, C=C)
    
def beta():
    return Expression(('1.0','0.0'))

def u0_boundary(x, on_boundary):
	return on_boundary
def outflow(x):
	return abs(x[0]-1) < 1E-14 or abs(x[1]-1) < 1E-14 or abs(x[1]) < 1E-14
def inflow(x):
	return abs(x[0]) < 1E-14 #or abs(x[1]-1) < 1E-14 or abs(x[1]) < 1E-14
class Inflow(Expression):
	def eval(self, values, x):
		values[0] = 0.0
		if inflow(x):
			values[0] = 1.0
