import getopt, sys

# ================= dump to matlab ========================= 

def dumpMat(name, AA):
    f = open(name, 'w')
    f.write("%s = sparse(%d,%d);\n" % (name,AA.shape[0],AA.shape[1]))
    for i in range(AA.shape[0]):
        for j in range(AA.shape[1]):
           if abs(AA[i,j]) > 10e-10:
               f.write("%s (%d, %d) = %e;\n " % (name,i+1,j+1,AA[i,j]))

def dump_i_vec(name, b):
	f = open(name, 'w')
	for i in range(len(b)):
		f.write("%s (%d) = %d;\n " % (name,i+1,b[i]))

# double precision vector
def dump_d_vec(name, b):
	f = open(name, 'w')
	for i in range(len(b)):
		f.write("%s (%d) = %e;\n " % (name,i+1,b[i]))

# example: parseArg('--eps',sys,default_val)                
def parseArg(*arg):
    argName = arg[0]
    sys = arg[1]
    default_val = ''
    if len(arg)>2:
        default_val = arg[2]

    argvec = ['eps=','N=','p=','numRefs=']
    inputs, remainder = getopt.getopt(sys.argv[1:],'',argvec)
    #print sys.argv
    #print inputs
    inputFound = False
    for arg,val in inputs:
        if arg==argName:
            inputFound=True
            return val
    if (inputFound==False):
        return default_val

eps = float(parseArg('--eps',sys,1e-2))
print eps
