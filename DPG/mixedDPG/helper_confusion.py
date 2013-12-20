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
    return Expression('(exp(r1*(x[0]-1))-exp(r2*(x[0]-1)))/(exp(-r1)-exp(-r2))*sin(pi*x[1])',eps = eps,pi = 2*acos(0.0), r1 = (-1+sqrt(1+4*pi*pi*eps*eps))/(-2*eps),r2 = (-1-sqrt(1+4*pi*pi*eps*eps))/(-2*eps))
    
class Inflow(Expression):
	def eval(self, values, x):
		values[0] = 0.0
		if abs(x[0]) < 1E-14:
			values[0] = 1.0
