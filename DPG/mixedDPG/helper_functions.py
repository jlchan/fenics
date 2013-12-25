import getopt, sys

# ================= dump to matlab ========================= 

def dumpMat(name, AA):
    f = open("read"+name+".m", 'w')
    f.write("%s = sparse(%d,%d);\n" % (name,AA.shape[0],AA.shape[1]))
    for i in range(AA.shape[0]):
        for j in range(AA.shape[1]):
           if abs(AA[i,j]) > 10e-10:
               f.write("%s (%d, %d) = %e;\n " % (name,i+1,j+1,AA[i,j]))

def dump_i_vec(name, b):
	f = open("read"+name+".m", 'w')
	for i in range(len(b)):
		f.write("%s (%d) = %d;\n " % (name,i+1,b[i]))

# double precision vector
def dump_d_vec(name, b):
	f = open(name, 'w')
	for i in range(len(b)):
		f.write("%s (%d) = %e;\n " % (name,i+1,b[i]))

# example: parseArg('--eps',default_val)
def parseArg(*arg):
    argName = arg[0]
    default_val = ''
    if len(arg)>1:
        default_val = arg[1]
        
    # add all arguments 
    argvec = [a.split("--")[1]+"=" for a in sys.argv[1:] if '--' in a]
    inputs, remainder = getopt.getopt(sys.argv[1:],'',argvec)

    inputFound = False
    for arg,val in inputs:
        if arg==argName:
            inputFound = True
            return val
    return default_val

