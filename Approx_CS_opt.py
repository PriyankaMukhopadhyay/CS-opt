import numpy as np
import math
from math import pi
import time
import random

n=2
N=2**n
N2=4**n
epsilon = 0.01
r_dig = 2 #Rounding place

# Define the matrices
I = np.array([[1, 0], [0, 1]])
X = np.array([[0, 1], [1, 0]])
Y = np.array([[0, -1j], [1j, 0]])
Z = np.array([[1, 0], [0, -1]])
T = np.array([[1,0],[0, np.exp((1j*pi)/4)]])

def tensor_product_recursive(paulis, depth):
        if depth == 1:
            return paulis
        else:
            return [np.kron(p, q) for p in paulis for q in tensor_product_recursive(paulis, depth - 1)]

def tensor_product(operators):
    result = operators[0]
    for op in operators[1:]:
        result = np.kron(result, op)
    return result

def trace_product(W_prime, P_out, P_in):
    product = np.matmul(W_prime, np.matmul(P_out, np.matmul(W_prime.conj().T, P_in)))
    return np.trace(product)

def matrix_product_recursive(S, depth):
        if depth == 1:
            return S
        else:
            return [np.matmul(p, q) for p in S for q in matrix_product_recursive(S, depth - 1)]

def matrix_neg_check(A, B, N):
    for i in range(N):
        for j in range(N):
            if A[i][j] != -B[i][j]:
                return False
    return True

def matrix_eq_check(A,B,N): #Returns 1 if equal
    for i in range(N):
      if eq == 0:
        break
      for j in range(N):
        #if (A[i][j] - B[i][j] < 10**-6) or (B[i][j] - A[i][j] < 10**-6):
        if (A[i][j] - B[i][j] == 0):
          eq = 1
        else:
          eq = 0
          break
    return eq

def generate_pauli_n(n):
    """Generate the set of Pauli matrices for n qubits."""
    # Base Pauli matrices
    I = np.array([[1, 0], [0, 1]])
    X = np.array([[0, 1], [1, 0]])
    Y = np.array([[0, -1j], [1j, 0]])
    Z = np.array([[1, 0], [0, -1]])
    T = np.array([[1,0],[0, np.exp((1j*pi)/4)]])
    paulis = [I, X, Y, Z]

    return tensor_product_recursive(paulis, n)

def pauli_string_to_matrix(pauli_string):
    # Convert a Pauli string (e.g., 'IXY') to its matrix representation
    # You might need to implement this function based on your requirements

    pauli_dict = {'I': I, 'X': X, 'Y': Y, 'Z': Z}
    matrix = np.kron(pauli_dict[pauli_string[0]], pauli_dict[pauli_string[1]])

    for p in pauli_string[2:]:
      matrix = np.kron(matrix, pauli_dict[pauli_string[2]])

    return matrix

pauli_I = [I]
i_tensored = tensor_product_recursive(pauli_I, n)[0]
pauli_n = generate_pauli_n(n)

def CLIFF_TEST(W_prime): #Returns 1 if Clifford

  i = 0
  cliff = 1
  for P in pauli_n:
    result = np.matmul(W_prime, P)
    tr_result = np.abs(np.trace(result) / N)
    if (tr_result != 0) and (i == 0):
      if tr_result == 1:
        break
      else:
        prev_amp = tr_result
        i = i+1
    if (tr_result != 0) and (i != 0):
      if tr_result != prev_amp:
        print("Not Clifford after amplitude test.")
        cliff = 0
        return 0

  for P_out in pauli_n:
    if cliff == 0:
      print("Exiting outer loop of conjugation test")
      return 0
    if np.allclose(P_out, i_tensored) == True:
      #print("Got I!")
      continue

    for P_in in pauli_n:
      if np.allclose(P_in, i_tensored) == True:
        #print("Got I!")
        continue
      val = abs(trace_product(W_prime, P_out, P_in)) / N
      if val == 1:
        break
      if (val != 1) and (val != 0):
        print("Not Clifford after conjugation test.")
        cliff = 0
        break

  if cliff == 1:
    print("Passed Clifford test")
    return 1

b_conj = 1 - 4 * epsilon**2 + 2 * epsilon**4
print("b_conj = ",b_conj)
#b_conj = round(b_conj,r_dig)
print("b_conj after rounding  = ",b_conj)

l_conj = 2 * epsilon
print("l_conj = ",l_conj)
#l_conj = round(l_conj,r_dig)
print("l_conj after rounding = ",l_conj)

def A_CONJ(W_prime, epsilon):

    p = 1

    for P_out in pauli_n:
        if np.allclose(P_out, i_tensored) == True:
            #print("Got I!")
            continue
        if p == 1:
            p = 0

        prod_out = np.matmul(W_prime, np.matmul(P_out, W_prime.conj().T))
        for P_in in pauli_n:
            if np.allclose(P_in, i_tensored) == True:
                #print("Got I!")
                continue
            prod_in = np.matmul(prod_out, P_in)
            val0 = abs(np.trace(prod_in)) / N
            val = round(val0,r_dig)
            print("val0,val, b_conj,l_conj = ",val0,val, b_conj,l_conj)
            #b_conj = 1 - 4 * epsilon**2 + 2 * epsilon**4

            if b_conj <= val <= 1:
                p += 1
                print("satisfies 1st ineq: p = ",p)
                if p > 1:
                    return "NO"
            else:
                print("notsatisfies 1st ineq")

            if l_conj < val < b_conj:
                print("satisfies 2nd ineq")
                return "NO"

    return "YES"

LB1 = []
UB1 = []
UB0 = []

for M in range(1,N2+1):
  lb_1 = (1 - epsilon**2) / np.sqrt(M) - np.sqrt(M * (2 * epsilon**2 - epsilon**4))
  LB1.append(lb_1)
  ub_1 = (1 / np.sqrt(M)) + np.sqrt(M * (2 * epsilon**2 - epsilon**4))
  UB1.append(ub_1)
  ub_0 = np.sqrt(M * (2 * epsilon**2 - epsilon**4))
  UB0.append(ub_0)

print("LB1 = ",LB1)
print("UB1 = ",UB1)
print("UB0 = ",UB0)

#for M in range(0,N2):
  #LB1[M] = round(LB1[M],r_dig)
  #UB1[M] = round(UB1[M],r_dig)
  #UB0[M] = round(UB0[M],r_dig)

print("LB1 after rounding = ",LB1)
print("UB1 after rounding = ",UB1)
print("UB0 after rounding = ",UB0)

def A_DECIDE(W,U_tilde, epsilon):

  W_prime = np.matmul(np.conj(W.T), U_tilde)
  #tr_Wprime = np.abs(np.trace(W_prime) / N)
  #dist = math.sqrt(1-tr_Wprime)
  S_c = []

  for P in pauli_n:
    result = np.matmul(W_prime, P)
    tr_result = np.abs(np.trace(result) / N)
    tr_result = round(tr_result,r_dig)
    S_c.append(tr_result)

  S_c.sort(reverse=True)

  for M in range(1,4**n+1):
      S_1 = S_c[:M]
      S_0 = S_c[M:]
      #lb1 = (1 - epsilon**2) / np.sqrt(M) - np.sqrt(M * (2 * epsilon**2 - epsilon**4))
      #ub1 = (1 / np.sqrt(M)) + np.sqrt(M * (2 * epsilon**2 - epsilon**4))
      #ub0 = np.sqrt(M * (2 * epsilon**2 - epsilon**4))
      amp = 1
      for x in S_1:
        if LB1[M-1] <= x <= UB1[M-1]:
          amp = 1
        else:
          amp = 0
          break
        if amp == 1:
          for x in S_0:
            if 0 <= x <= UB0[M-1]:
              amp = 1
            else:
              amp = 0
              break
        if amp == 1:
          print("Conj begin")
          result = A_CONJ(W_prime, epsilon)
          if result == "YES":
            execTime = time.time()-start_time
            #print("Fin Utilde",U_tilde)
            #print("Fin S_c",S_c)
            #print("Fin M",M)
            #print("Fin S_1",S_1)
            #print("Fin S_0",S_0)
            return 2 #Passed both tests
          else:
            return 1 # Passed amp, did not pass conj

  return 0

#-----------------------SETS-----------------------------
S1 = []
S2 = []
S3 = []
S4 = []
S5 = []
S6 = []
S7 = []
S8 = []
S9 = []
S10 = []

cs_gen = [('IX', 'XI'), ('IX', 'YI'), ('IX', 'ZI'), ('IY', 'XI'), ('IY', 'YI'), ('IY', 'ZI'), ('IZ', 'XI'), ('IZ', 'YI'), ('IZ', 'ZI'), ('XX', 'YY'), ('XX', 'YZ'), ('XY', 'YX'), ('XY', 'YZ'), ('XZ', 'YX'), ('XZ', 'YY')]
cs_gen_t = []
beta_0 = (3 + 1j)/4
beta_1 = (1 - 1j)/4
for l in range(len(cs_gen)):
  P1 = pauli_string_to_matrix(cs_gen[l][0])
  P2 = pauli_string_to_matrix(cs_gen[l][1])
  G_P12 = beta_0*i_tensored + beta_1*(P1 + P2 - np.matmul(P1, P2))
  cs_gen_t.append(G_P12)
  #print("l,P1,P2,cs_gen_t = ",l,P1,P2,cs_gen_t[l])

print("size of cs_gen_t = ",len(cs_gen_t))
size_gcs = len(cs_gen)
print("size_gcs = ",size_gcs)

for l1 in range(size_gcs):
  prod1 = cs_gen_t[l1]
  S1.append(prod1)
  for l2 in range(size_gcs):
    if l2 == l1:
      continue
    prod2 = np.matmul(prod1,cs_gen_t[l2])
    S2.append(prod2)
    for l3 in range(size_gcs):
      if l3 == l2:
        continue
      prod3 = np.matmul(prod2,cs_gen_t[l3])
      S3.append(prod3)
      for l4 in range(size_gcs):
        if l4 == l3:
          continue
        prod4 = np.matmul(prod3,cs_gen_t[l4])
        S4.append(prod4)
        for l5 in range(size_gcs):
          if l5 == l4:
            continue
          prod5 = np.matmul(prod4,cs_gen_t[l5])
          S5.append(prod5)
          for l6 in range(size_gcs):
            if l6 == l5:
              continue
            prod6 = np.matmul(prod5,cs_gen_t[l6])
            S6.append(prod6)
            #for l7 in range(size_gcs):
              #if l7 == l6:
                #continue
              #prod = np.matmul(prod,cs_gen_t[l7])
              #S7.append(prod)

size_S1 = len(S1)
size_S2 = len(S2)
size_S3 = len(S3)
size_S4 = len(S4)
size_S5 = len(S5)
size_S6 = len(S6)
#size_l7 = len(S7)

print("Size of S1 = ",size_S1)
print("Size of S2 = ",size_S2)
print("Size of S3 = ",size_S3)
print("Size of S4 = ",size_S4)
print("Size of S5 = ",size_S5)
print("Size of S6 = ",size_S6)
#print("Size of S7 = ",size_S7)

#--------------------------------UNITARY-----------------------
#W = np.array([[1,0],[0, np.exp((1j*pi)/4)]])  #T gate
theta_k = pow(2,2)  #k=0 : Z, k=1: S, k=2 : T
W = np.array([[1,0,0,0],[0, np.exp((1j*pi)/theta_k),0,0],[0,0,1,0],[0,0,0,np.exp((1j*pi)/theta_k)]]) # Example 2x2 unitary matrix (identity matrix)
print("W, epsilon = ",W, epsilon)
#epsilon = 0.05

cliff_result = CLIFF_TEST(W)
if cliff_result == 1:
  print("Found Clifford")
else:
  print("Not Clifford")

#---------------------OPT--------------------
start_time = time.time()
succ = 0
i1 = 0

while i1  < size_S6:
  i1 = i1+1
  if succ == 1:
    break
  U_1 = S6[i1]
  for i2 in range(size_S6):
    U_tilde = np.matmul(U_1, S6[i2])
    result = A_DECIDE(W,U_tilde, epsilon)
    if result == "YES":
      print("YES")
      succ = succ+1
      break


execTime = time.time()-start_time
print("Implementation time = ",execTime)




