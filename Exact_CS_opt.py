import numpy as np
import time
import random
import math

n=3
N=2**n
N2=4**n

# Define the matrices
I = np.array([[1, 0], [0, 1]], dtype=complex) #intg = 0
X = np.array([[0, 1], [1, 0]], dtype=complex) #intg = 1
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)  #intg = 2
Z = np.array([[1, 0], [0, -1]], dtype=complex)  #intg = 3

pauli_1=[I,X,Y,Z]

def MULTIQUBIT_PAULI(qubit):
    if qubit == 1:
      return pauli_1
    else:
      P_n = []
      for p in pauli_1:
        for q in MULTIQUBIT_PAULI(qubit-1):
          temp=np.kron(p,q)
          P_n.append(temp)
    return P_n

pauli_n = MULTIQUBIT_PAULI(n)
print("length of pauli_n = ",len(pauli_n))

def quartDecompose(N_func):
    N_decompose=[0]*n
    i = n-1
    #m_temp=1
    m=4
    for j in range(n):
      #m=pow(4,m_temp)
      N_decompose[i]=int(N_func%m)
      N_func=(N_func-N_decompose[i])/4
      i=i-1
    return N_decompose

def quartInt(N_arr):  #Arr is indexed as 0123...(n-1)
    sum = 0
    for j in range(n):
      sum = sum + N_arr[n-1-j]*pow(4,j)

    return sum

def ARR_MATRIX(P):
    for j in range(n):
      if j == 0:
        mat = pauli_1[P[0]]
      else:
        mat=np.kron(pauli_1[P[j]],mat)

    return mat

def ARR_EQ(A,B,size):  #Returns 0 if equal
    eq = 0
    for i in range(size):
      if A[i] != B[i]:
        eq = 1
        break

    return eq

def PAULI_MULT(P, Q):  #No phase
    R = [0] * n  # Initialize the result with identity Pauli
    #r = 1  # Initialize r to +1

    for j in range(n):
        if P[j] == 0:
            R[j] = Q[j]
        elif Q[j] == 0:
            R[j] = P[j]
        elif P[j] == Q[j]:
            R[j] = 0
        elif P[j] == 1 and Q[j] == 2:
            R[j] = 3
            #r *= 1j
        elif P[j] == 2 and Q[j] == 1:
            R[j] = 3
            #r *= -1j
        elif P[j] == 2 and Q[j] == 3:
            R[j] = 1
            #r *= 1j
        elif P[j] == 3 and Q[j] == 2:
            R[j] = 1
            #r *= -1j
        elif P[j] == 3 and Q[j] == 1:
            R[j] = 2
            #r *= 1j
        elif P[j] == 1 and Q[j] == 3:
            R[j] = 2
            #r *= -1j

    #R.insert(0, r)
    return R

def PAULI_COMM(P, Q):
    x = 0

    for j in range(n):
        if (P[j] != 0) and (Q[j] != 0) and (P[j] != Q[j]):
            x = x+1

    if x % 2 == 0:
        return 1
    else:
        return 0

#---------------calculate G_CS as array of pairs of Paulis-----------------

G_CS = []
for P1_int in range(1,N2):
  P1 = quartDecompose(P1_int)
  for P2_int in range(1,N2):
    if P2_int <= P1_int:
      continue
    P2 = quartDecompose(P2_int)
    if PAULI_COMM(P1,P2) == 1 :
      flag = 1
      #print("P1,P2 = ",P1,P2)
      for [Pa_int,Pb_int] in G_CS:
        Pa = quartDecompose(Pa_int)
        Pb = quartDecompose(Pb_int)
        prod = PAULI_MULT(Pa,Pb)
        prod_int = quartInt(prod)
        flag_match = 0
        R = [Pa_int,Pb_int]
        R_prime = [P1_int,P2_int]
        Q=[0,0]
        indx1=[0,1]
        indx2=[0,1]
        for j in range(2):
          for k in range(2):
            if R[j] == R_prime[k]:
              Q[j] = R_prime[k]
              flag_match = 1
              indx1.remove(j)
              indx2.remove(k)
        if flag_match == 1:
          Q[indx1[0]] = R_prime[indx2[0]]
          if (prod_int == Q[0]) or (prod_int == Q[1]):
            flag = 0
            break
      if flag == 1:
        G_CS.append([P1_int,P2_int])

print("G_CS = ",G_CS)
print("Size of G_CS = ",len(G_CS))

#-------------CHANNEL REP AS ARRAY---------------------------

chanG_CS = []

for [P1_int,P2_int] in G_CS:
  A_P1P2 = []
  P1 = quartDecompose(P1_int)
  P2 = quartDecompose(P2_int)
  P1_mat = ARR_MATRIX(P1)
  P2_mat = ARR_MATRIX(P2)
  Q = P1_mat + P2_mat - np.matmul(P1_mat,P2_mat)
  for r in range(N2):
    r_decompose = quartDecompose(r)
    Pr = ARR_MATRIX(r_decompose)
    #Pr = pauli_n[r]
    tup = []
    mat_temp = np.matmul(np.matmul(Pr,Q),np.matmul(Pr,Q))
    val =(1/N)*( (5/8)*np.trace(pauli_n[0])+(1/8)*np.trace(mat_temp) )
    val_real = val.real
    #print("val_real = ",val_real)
    if val_real == 0.5:
      tup.append(r)
      num = 0
      for s in range(N2):
        if s == r:
          continue
        s_decompose = quartDecompose(s)
        Ps = ARR_MATRIX(s_decompose)
        #Ps = pauli_n[s]
        val =  (1/N)*( ((1+2j)/8)*np.trace(np.matmul(np.matmul(Pr,Ps),Q )) + ((1-2j)/8)*np.trace(np.matmul(np.matmul(Ps,Pr),Q ))
        + (1/8)*np.trace( np.matmul(np.matmul(Pr,Q),np.matmul(Ps,Q)) ) )
        val_real = val.real
        #print("val_real = ",val_real)
        if val_real == 0.5:
          tup.append(s)
          num = num+1
          if num == 3:
            A_P1P2.append(tup)
            break
        if val_real == -0.5:
          tup.append(-s)
          num = num+1
          if num == 3:
            A_P1P2.append(tup)
            break

  #print("[P1,P2] = ",P1_int,P2_int)
  #print("A_P1P2 = ",A_P1P2)
  #print("Size of A_P1P2=",len(A_P1P2))
  chanG_CS.append(A_P1P2)

print("chanG_CS = ",chanG_CS)

size_GCS = len(G_CS)

print("size of G_CS = ",size_GCS)
print("size of chanG_CS = ",len(chanG_CS))

#-------------CHANNEL REP INVERSE AS ARRAY---------------------------

chanG_CSinv = []

for [P1_int,P2_int] in G_CS:
  A_P1P2 = []
  P1 = quartDecompose(P1_int)
  P2 = quartDecompose(P2_int)
  P1_mat = ARR_MATRIX(P1)
  P2_mat = ARR_MATRIX(P2)
  Q = P1_mat + P2_mat - np.matmul(P1_mat,P2_mat)
  for r in range(N2):
    Pr = pauli_n[r]
    tup = []
    mat_temp = np.matmul(np.matmul(Pr,Q),np.matmul(Pr,Q))
    val =(1/N)*( (5/8)*np.trace(pauli_n[0])+(1/8)*np.trace(mat_temp) )
    val_real = val.real
    #print("val_real = ",val_real)
    if val_real == 0.5:
      tup.append(r)
      num = 0
      for s in range(N2):
        if s == r:
          continue
        Ps = pauli_n[s]
        val =  (1/N)*( ((1-2j)/8)*np.trace(np.matmul(np.matmul(Pr,Ps),Q )) + ((1+2j)/8)*np.trace(np.matmul(np.matmul(Ps,Pr),Q ))
        + (1/8)*np.trace( np.matmul(np.matmul(Pr,Q),np.matmul(Ps,Q)) ) )
        val_real = val.real
        #print("val_real = ",val_real)
        if val_real == 0.5:
          tup.append(s)
          num = num+1
          if num == 3:
            A_P1P2.append(tup)
            break
        if val_real == -0.5:
          tup.append(-s)
          num = num+1
          if num == 3:
            A_P1P2.append(tup)
            break

  #print("[P1,P2] = ",P1_int,P2_int)
  #print("A_P1P2 = ",A_P1P2)
  #print("Size of A_P1P2=",len(A_P1P2))
  chanG_CSinv.append(A_P1P2)

print("chanG_CSinv = ",chanG_CSinv)

print("size of G_CS = ",size_GCS)
print("size of chanG_CSinv = ",len(chanG_CSinv))


#------------SUB-ROUTINES FOR EXACT-CS-DECIDE-----------

def GET_SDE(U):

  sde=0

  for i in range(N2):
    for j in range(N2):
      if U[i][j][1] > sde:
          sde = U[i][j][1]

  return sde
  
  def sde2_REDUCE(a, k):
    while a % 2 == 0 and a != 0 and k!=0:  # Check if 'a' is divisible by 2 and not zero
        a = a // 2  # Use integer division to avoid floating point result
        #print(a)
        k = k - 1
        #print(k)
    return [a,k]

def ADD_2(v1, v2):
    [a1, k1], [a2, k2] = v1, v2

    if k1 >= k2:
        num = a1 + a2 * (2**(k1 - k2))
        if num == 0:
          den = 0
        else:
          den = k1
    else:
        num = a1 * (2**(k2 - k1)) + a2
        if num == 0:
          den = 0
        else:
          den = k2

    return sde2_REDUCE(num, den)

def MULT_GCS(Pindx, U, inv): #inv=1 if inverse
    # Input P is a number from 0 to 4^n. s_in is sde of input matrix U.
    Up = [[[0, 0] for _ in range(N2)] for _ in range(N2)]

    for i in range(N2):
      for j in range(N2):
        Up[i][j] = U[i][j]

    if inv == 0:
      A_P1P2 = chanG_CS[Pindx]
    else:
      A_P1P2 = chanG_CSinv[Pindx]


    for i in range(int((3*N2)/4)):
      diag = A_P1P2[i][0]
      oDiag1 = A_P1P2[i][1]
      oDiag2 = A_P1P2[i][2]
      oDiag3 = A_P1P2[i][3]
      
      if oDiag1 < 0:
        oDiagIndx1 = -oDiag1
      else:
        oDiagIndx1 = oDiag1
      if oDiag2 < 0:
        oDiagIndx2 = -oDiag2
      else:
        oDiagIndx2 = oDiag2
      if oDiag3 < 0:
        oDiagIndx3 = -oDiag3
      else:
        oDiagIndx3 = oDiag3
      for j in range(N2):
        v1 = U[diag][j]
        v2 = U[oDiagIndx1][j]
        v3 = U[oDiagIndx2][j]
        v4 = U[oDiagIndx3][j]
        
        if v1[0] != 0:
          v1 = [v1[0], v1[1]+1]
        if v2[0] != 0:
          if oDiag1 < 0:
            v2 = [-v2[0], v2[1]+1]
          else:
            v2 = [v2[0], v2[1]+1]
        if v3[0] != 0:
          if oDiag2 < 0:
            v3 = [-v3[0], v3[1]+1]
          else:
            v3 = [v3[0], v3[1]+1]
        if v4[0] != 0:
          if oDiag3 < 0:
            v4 = [-v4[0],v4[1]+1]
          else:
            v4 = [v4[0], v4[1]+1]
        
        sum12=ADD_2(v1,v2)
        sum34=ADD_2(v3,v4)
        
        Up[diag][j] = ADD_2(sum12,sum34)
      

    return Up


def HAM_WT_MAT(U):
    ham = 0
    for i in range(N2):
        for j in range(N2):
            if U[i][j][0] != 0:
                ham += 1
    return ham
    
def UPDATE_SH(SH, sde1, ham1, sde0, ham0, rule):
    #s = h = 0  # Initialize s and h to 0 sde1 : child, sde0 : par
    # inc (0), unchanged (1), dec (2) : sde along row and ham along column
    #SH is always 3x3

    if rule == 1: # 9 groups : Both sde and ham inc, unchanged, dec
      if sde1 > sde0 and ham1 > ham0:
        SH[0][0] = SH[0][0] + 1; s = h = 0
      elif sde1 > sde0 and ham1 == ham0:
        SH[0][1] = SH[0][1] + 1; s = 0; h = 1
      elif sde1 > sde0 and ham1 < ham0:
        SH[0][2] = SH[0][2] + 1; s = 0; h = 2
      elif sde1 == sde0 and ham1 > ham0:
        SH[1][0] = SH[1][0] + 1; s = 1; h = 0
      elif sde1 == sde0 and ham1 == ham0:
        SH[1][1] = SH[1][1] + 1; s = 1; h = 1
      elif sde1 == sde0 and ham1 < ham0:
        SH[1][2] = SH[1][2] + 1; s = 1; h = 2
      elif sde1 < sde0 and ham1 > ham0:
        SH[2][0] = SH[2][0] + 1; s = 2; h = 0
      elif sde1 < sde0 and ham1 == ham0:
        SH[2][1] = SH[2][1] + 1; s = 2; h = 1
      elif sde1 < sde0 and ham1 < ham0:
        SH[2][2] = SH[2][2] + 1; s = 2; h = 2

    if rule == 4: #2 groups : sde inc, non-inc (same+dec)
      if sde1 > sde0:
        SH[0][0] = SH[0][0] + 1; s = h = 0
      elif sde1 <= sde0:
        SH[2][0] = SH[2][0] + 1; s = 2; h = 0

    return (SH, s, h)
    
def MIN_SH(SH):

  start = 0
  for i in range(3):
    for j in range(3):
      if (start == 0) and (SH[i][j] > 0) :
        start = 1
        min_val = SH[i][j]
        s_index = i
        h_index = j
      if (start == 1) and (SH[i][j] > 0):
        if SH[i][j] <= min_val:
          min_val = SH[i][j]
          s_index = i
          h_index = j

  return (s_index, h_index)


#---------------EXACT-CS-DECIDE------------------

def EXACT_CS_DECIDE(U, m, sde_U):

  ham_U=HAM_WT_MAT(U)
  Path_U=[]
  Par_node=[]
  U_tilde=[U,Path_U,sde_U,ham_U]
  #print("U_tilde",U_tilde)
  Par_node.append(U_tilde)
  print("Root = ",Par_node)
  SH=np.zeros((3,3),dtype=int)
  rule = 1
  leaf = 0
  max_sel_node = 1

  for i in range(1,m+1):
    if leaf == 1:
      #print("Breaking due to leaf==1 at i = ",i)
      break
    Child_node=[]
    print("Level = ",i)
    num_par = len(Par_node)
    if num_par > max_sel_node:
      max_sel_node = num_par
    #print(" No of Parents = ",len(Par_node))
    SH=np.zeros((3,3),dtype=int)
    if num_par == 0:
      #print("No of parent nodes is ..so breaking",len(Par_node))
      break
    for j in range(num_par):
      if leaf == 1:
        #print("Breaking due to leaf ==1 at j = ",j)
        break
      U_par=Par_node[j][0]
      #print("U+par = ",U_par)
      Path_par=Par_node[j][1]
      print("Path_par = ",Path_par)
      sde_par=Par_node[j][2]
      ham_par=Par_node[j][3]
      path_len=len(Path_par)
      if path_len > 0:
        P_prev=Path_par[path_len-1]
      for k in range(size_GCS):
        Path_W=[]
        P=k
        if (path_len > 0) and (P == P_prev):
          continue
        if path_len > 0:
          for p_c in range(path_len):
            Path_W.append(Path_par[p_c])
        #print("P=",P)
        W=MULT_GCS(P,U_par,1)
        sde_W=GET_SDE(W)
        print("P,sde_W = ",P,sde_W)
        ham_W=HAM_WT_MAT(W)
        Path_W.append(P)
        if sde_W == 0:
          fin_i = i
          fin_Path = Path_W
          leaf = 1
          break
        if sde_W == 1:
          W_temp =[[[0,0] for i in list(range(N2))] for j in list(range(N2))]
          for P_temp in range(size_GCS):
            W_temp = MULT_GCS(P_temp,W,1)
            sde_Wtemp = GET_SDE(W_temp)
            if sde_Wtemp == 0:
              fin_i = i+1
              Path_W.append(P_temp)
              fin_Path = Path_W
              leaf = 1
              print("Breaking at sde 1")
              break
            else:
              print("No breaking at sde 1")
          if leaf == 1:
            print("Exiting outer for loop")
            break
        #print("P,sde_W,m+1-i, case 2=",P,sde_W,m+1-i)
        #continue
        if (sde_W > 0) and (leaf == 0):
          new_SH=UPDATE_SH(SH,sde_W,ham_W,sde_par,ham_par,rule)
          #print("P,sde_W,sde_par,m+1-i, Case 3 = ",P,sde_W,sde_par,m+1-i)
          SH=new_SH[0]
          s_W=new_SH[1]
          h_W=new_SH[2]
          Child_node.append([W,s_W,h_W,Path_W,sde_W,ham_W])

    num_child=len(Child_node)
    if (leaf == 0) and (num_child != 0):
      #print("Children = ",Child_node)
      #print("No of children = ",len(Child_node))
      sh=MIN_SH(SH)
      s_indx=sh[0]
      h_indx=sh[1]
      #print("SH,s_indx,h_index = ",SH,s_indx,h_indx)
    Par_node = []

    if num_child != 0:
      for j in range(len(Child_node)):
        if leaf == 1:
          print("breaking due to leaf == 1 at child j,",j)
          break
        if Child_node[j][4] > m+1-i:
            #print("P,sde_W,m+1-i, case 2=",P,sde_W,m+1-i)
            continue
        if Child_node[j][4] == 1:
          next_U=Child_node[j][0]
          next_U_Path=Child_node[j][3]
          next_U_sde=Child_node[j][4]
          next_U_ham=Child_node[j][5]
          Par_node.append([next_U,next_U_Path,next_U_sde,next_U_ham])
        if (Child_node[j][1] == s_indx) and (Child_node[j][2] == h_indx):
          next_U=Child_node[j][0]
          next_U_Path=Child_node[j][3]
          next_U_sde=Child_node[j][4]
          next_U_ham=Child_node[j][5]
          Par_node.append([next_U,next_U_Path,next_U_sde,next_U_ham])

  if leaf == 1:
    print("Max no of selected nodes = ",max_sel_node)
    return (fin_i,fin_Path)
  else :
    return (-1,[])



#-------------CHAN REP FOR Toffoli--------------------
#Code for specific unitaries ---To delete if implementing random unitaries, defined next.

#Code for specific unitaries ---

Tof =  np.array([[1, 0, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0, 1, 0]], dtype=complex)

print("Tof = ",Tof)
Tof_adj = Tof.conj().T
print("Tof_adj = ",Tof_adj)

U_in =[[[0,0] for i in list(range(N2))] for j in list(range(N2))]
diag = 0
noDiag = 0
poDiag = 0


for i in range(N2):
  Pr = pauli_n[i]
  #print("Pr = ",Pr)
  for j in range(N2):
    Ps = pauli_n[j]
    #print("Ps = ",Ps)
    prod = np.matmul(Pr,np.matmul(Tof,np.matmul(Ps,Tof_adj)))
    val = np.trace(prod)/N
    val_r = val.real
    if val_r == 1.0 :
      U_in[i][j] = [1,0]
    if val_r == 0.5:
      U_in[i][j] = [1,1]
      if i == j:
        diag = diag+1
      else:
        poDiag = poDiag+1
    if val_r == -0.5:
      U_in[i][j] = [-1,1]
      if i == j:
        print("Err, i, j, val_r = ",i,j,val_r)
      else:
        noDiag = noDiag+1

print("U_in = ",U_in)
print("diag, poDiag, noDiag, sum = ",diag, poDiag, noDiag, diag+poDiag+noDiag)




#------------------RANDOM-CHAN-REP---------------
# Code for Random Unitaries. To delete if implementing specific unitaries

# Code for Random Unitaries
cs_in = 10

U_temp =[[[0,0] for i in list(range(N2))] for j in list(range(N2))]
for i in range(N2):
  U_temp[i][i]=[1,0]

print("U_temp = ",U_temp)
sde_in=0
i=0

while i < cs_in:
  P_in = np.random.randint(0, size_GCS)
  print("P_in = ",P_in)

  if i == 0:
    P_prev = P_in
    U_temp = MULT_GCS(P_in, U_temp,0)
    i = i+1
  else:
    if P_in == P_prev:
      continue
    else:
      P_prev = P_in
      U_temp = MULT_GCS(P_in, U_temp,0)
      i = i+1

#U_temp = MULT_GCS(P_in, U_temp, 1)

sde_in=GET_SDE(U_temp)
print("Final sde = ",sde_in)
print("U_temp = ",U_temp)

cols = list(range(0,N2))
print("cols = ",cols)
random.shuffle(cols)
print("After perm cols = ",cols)

U_in =[[[0,0] for i in list(range(N2))] for j in list(range(N2))]

for j in range(0,N2):
  temp_col = cols[j]
  b_in = np.random.randint(0, 2)
  #print("temp_col,j,b_in = ",temp_col,j,b_in)
  for i in range(0,N2):
    U_in[i][j] = U_temp[i][temp_col]
    if b_in == 1:
      U_in[i][j][0] = -U_in[i][j][0]

print("U_in = ",U_in)


#-------------------------EXACT-CS-OPT

#EXACT-CS-OPT
  #Input channel rep
print("Input Unitary = ", U_in)
sde_in = GET_SDE(U_in)
print("Input sde = ", sde_in)
m = sde_in
reach = 0

if m == 0:
    print("Clifford")

start_time = time.time()
while reach == 0:
    m_prime, D = EXACT_CS_DECIDE(U_in,m,sde_in)
    print("m_prime, D = ",m_prime,D)
    if m_prime == -1:
      m = m+1
      #reach = 1
    else:
      print("m_prime, D = ", m_prime,D)
      print("Path length = ",len(D))
      print("Input Unitary = ", U_in)
      reach = 1

execTime = time.time()-start_time
print("Time = ",execTime)

