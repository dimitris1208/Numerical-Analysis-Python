""" PA = LU factorization for root finding """

n = int(input("Type in the dimensions of the matrix A : "))
A = [ [0 for i in range(n) ] for j in range(n) ]
B = [  0 for i in range(n) ]
Y = [ 0 for i in range(n) ]
X = [ 0 for i in range(n) ]

print("Type the content of array A , row to row : ")
for i in range(n) :
    for j in range(n) :
        A[i][j] = float (input("A["+i+"]"+"["+j+"] <--- "))



print("Type the content of array B : ")
for i in range(n) :
    print("B["+i+"]")
    B[i] = float(input())



# The right row permutations happen in A and B matrix
for i in range(n) :
    if A[i][i] == 0 :
        temp1 = A[i].copy()
        temp2 = B[i]
        j=0
        done = False 
        while ( j < n and done == False) :
            if A[j][i] != 0 :
                A[i] = A[j]
                A[j] = temp1
                B[i] = B[j]
                B[j] = temp2
                done = True
            j += 1    



# PA = LU factorization with the changes being saved in the spots that A is zero
for i in range(n-1) :
    for j in range(i+1,n) :
        if A[j][i] != 0 :
            
            A[j][i] = ( A[j][i] / A[i][i] )
            for k in range(i+1,n) :
                A[j][k] = A[j][k] - ( A[j][i] * A[i][k] )


# LY = B is being solved , whereas L equals A's elements that are under the diagonal
for j in range(n) :
    Y[j] = B[j]
    for i in range(j) :
        Y[j] = Y[j] -A[j][i] * Y[i] 


# UX = Y is being solved , whereas U equals A's elements that are in the diagonal and above
for j in range(n-1,-1,-1) :
    X[j] = Y[j]
    for i in range(1,n-j) :
        X[j] = X[j]- A[j][j+i] * X[j+i]
    X[j] = X[j] / A[j][j]     

print(X)