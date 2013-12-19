from dolfin import *
import helper_functions as help
import helper_confusion as chelp

import sys,getopt

#print 'Number of arguments:', len(sys.argv), 'arguments.'
#print 'Argument List:', str(sys.argv)

p = int(help.parseArg('--p',sys,1))
N = int(help.parseArg('--N',sys,4))
numRefs = int(help.parseArg('--numRefs',sys,1))
eps = float(help.parseArg('--eps',sys,1e-2))

print "eps = ", eps

# define problem params
beta = Constant(('1.0','0.0'))
ue = chelp.erikkson_solution(eps)

def u0_boundary(x, on_boundary):
	return on_boundary

mesh = UnitSquareMesh(N,N)

# Load mesh and subdomains
h_vec = []
err_vec = []
for refIndex in xrange(numRefs):
    h = CellSize(mesh)
    
    # Create FunctionSpaces
    V = FunctionSpace(mesh, "CG", p)
    u,v = TrialFunction(V), TestFunction(V)
    
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
    tempMesh = chelp.quadrature_refine(mesh,N,4)    
    uh = interpolate(uh,FunctionSpace(tempMesh,"CG",p+2))
    e = sqrt(assemble((ue-uh)**2*dx))

    hh = (1.0/N)*.5**float(refIndex)
    h_vec.append(hh)
    err_vec.append(e)
    print "for h = ", hh, ", L2 error = ", e
    
    # refine mesh for next iteration
    mesh = refine(mesh)

plot(uh)
plot(tempMesh)
interactive()

help.dump_d_vec("h_eps"+str(eps)+"_p"+str(p)+".m",h_vec)
help.dump_d_vec("err_eps"+str(eps)+"_p"+str(p)+".m",err_vec)

