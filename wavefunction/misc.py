#
# I/O functions etc.
#
# J Robert Johansson, <robert@riken.jp>
#
from scipy import *
from scipy import linalg

#----
# Fig functions
#
def setFigure(Fig):

    if Fig:
        gplt.figure(Fig)
    else:
        gplt.figure()
        gplt.hold('on')
        Fig = gplt.current()
    
    return Fig


# -----------------------------------------------------------------------------------
#  Print functions           
#
#

#
# Print real part of matrix matrix to console
#
def print2DMatrix(A):
    N = A.shape[0]
    M = A.shape[1]

    for n in range(0, N):
        for m in range(0, M):
            val = real(A[n][m])
            if val > 0.0:
                print(" %.3f  " % (real(A[n][m]))),
            else:
                print("%.3f  " % (real(A[n][m]))),

        print



#
# Append row to file
#
def file_add_row(filename, row):

    file = open(filename, 'a')

    file.write('%s\n' % row)
    file.close

#
# Save a matrix in ASCII format to a file
#
def file_save_matrix(filename, M):

    n,m = shape(M)

    file = open(filename,'w')
    
    for i in range(0, n):
        for j in range(0, m):
            file.write('%.6f\t' % M[i,j])
        file.write('\n')
    
    file.close

#
# Save a row to a file
#
def file_save_row(filename, row):

    file = open(filename,'a')
    file.write(row + "\n")
    file.close

#
# Save vectors in ASCII format to a file
#
def file_save_vectors(filename, vectors, description, opt_real = False):

    vlen = len(vectors)
    mlen = 0

    file = open(filename,'w')
    
    if description:
        file.write('%s\n' % description)
    
    for v in arange(0, vlen):
        if mlen < len(vectors[v]):
            mlen = len(vectors[v])

    if opt_real:
        for i in arange(0, mlen):
            for v in arange(0, vlen):
                if len(vectors[v]) > i:
                    file.write('%.5E\t' % (real(vectors[v][i]))),
                else:
                    file.write('\t\t');
            file.write('\n')
        return
    

    for i in arange(0, mlen):
        for v in arange(0, vlen):
            if len(vectors[v]) > i:
                if (imag(vectors[v][i]) < 0):
                    file.write('%.5E%.5Ej\t' % (real(vectors[v][i]), imag(vectors[v][i]))),
                else: 
                    file.write('%.5E+%.5Ej\t' % (real(vectors[v][i]), imag(vectors[v][i]))),
            else:
                file.write('\t\t');# print('\t\t'),
        file.write('\n')

#
# Vector to CVS, TSV, etc.
#
def vector_2str(x, sep):
    rstr = ""
    for X in x:
        rstr += str(X)+sep
    return rstr
