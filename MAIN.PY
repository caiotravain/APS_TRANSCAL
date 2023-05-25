from funcoesTermosol import *
import numpy as np
import matplotlib.pyplot as plt

def matprint(mat, fmt="g"):
    col_maxes = [max([len(("{:"+fmt+"}").format(x)) for x in col]) for col in mat.T]
    for x in mat:
        for i, y in enumerate(x):
            print(("{:"+str(col_maxes[i])+fmt+"}").format(y), end="  ")
        print("")
  

def jacobi(ite, tol, K, F):
    tamanho = len(K)
    continuar = True

    erroMax = 0
    lista_x = []
    lista_err = []
    for i in range(tamanho):
        lista_x.append(0)
        lista_err.append(0)
        
    erroMax = 0

    numeroIteracoes = 0

    linhas_x = []
    for i in range(len(K)):
        linhas_x.append(0)
    
    while continuar:
        for i in range(len(linhas_x)):
            temp= 0
            temp2 = 0
            for j in range(len(linhas_x)):
                if i != j:
                    temp = K[i][j] * lista_x[j]
                    temp2 = temp + temp2            
            linhas_x[i] = (F[i][0] - temp2) / K[i][i]

       
            for i in range(len(linhas_x)):
                try:
                    if (linhas_x[i] != 0):
                        lista_err[i] = ((linhas_x[i] - lista_x[i]) / linhas_x[i])
                except:
                    pass
        

        for i in range(len(linhas_x)):
            lista_x[i] = linhas_x[i]

        numeroIteracoes += 1

        if numeroIteracoes == ite:
            continuar = False
        count = 0
        for i in range(len(lista_err)):
            if lista_err[i] < tol and numeroIteracoes > 1:
                count += 1
        if count == len(lista_err):
            continuar = False
     
    erroMax = max(lista_err)
    return [lista_x]



def gaussSeidel(ite, tol, K, F):
    continuar = True
    tamanho = len(K)

    erroMax = 0
    lista_x = []
    lista_err = []
    for i in range(tamanho):
        lista_x.append(0)
        lista_err.append(0)
        
    erroMax = 0

    numeroIteracoes = 0

    linhas_x = []
    for i in range(len(K)):
        linhas_x.append(0)
    
    while continuar:
        for i in range(len(linhas_x)):
            temp= 0
            temp2 = 0
            for j in range(len(linhas_x)):
                if i != j:
                    if i == 0:
                        temp = K[i][j] * lista_x[j]
                        temp2 = temp + temp2
                    else:
                        temp = K[i][j] * linhas_x[j]
                        temp2 = temp + temp2
                           
            linhas_x[i] = (F[i][0] - temp2) / K[i][i]

        try:
            for i in range(len(linhas_x)):
                if (linhas_x[i] != 0):
                    lista_err[i] = abs((linhas_x[i] - lista_x[i]) / linhas_x[i])
        except:
            pass
        

        for i in range(len(linhas_x)):
            lista_x[i] = linhas_x[i]

        numeroIteracoes += 1

        if numeroIteracoes == ite:
            continuar = False
        count = 0
        for i in range(len(lista_err)):
            if lista_err[i] < tol and numeroIteracoes > 1:
                count += 1
        if count == len(lista_err):
            continuar = False
     
    erroMax = max(lista_err)
    return [lista_x]


[nn,N,nm,Inc,nc,F,nr,R] = importa('entrada.xlsx')

plota(N,Inc)

lista_kg = np.zeros((nn*2, nn*2))

lista_deforma = []
lista_tensao = []
lista_ke = []



for i in range(nm):
    n1 = int(Inc[i,0])
    n2 = int(Inc[i,1])        
    E = Inc[i,2]
    A = Inc[i,3]
    x1 = N[0,n1-1]
    y1 = N[1,n1-1]
    x2 = N[0,n2-1]
    y2 = N[1,n2-1]
    L = np.sqrt((x2-x1)**2 + (y2-y1)**2)
    cos = (x2-x1)/L
    sen = (y2-y1)/L
    ke = np.array([[cos**2,cos*sen, -cos**2,-cos*sen], [cos*sen,sen**2,-cos*sen,-sen**2], [-cos**2,-cos*sen,cos**2,cos*sen], [-cos*sen,-sen**2,cos*sen,sen**2]])* E* A/L
    # print(ke)
    # print("---------------------------------")


    lista_ke.append(ke)


# Calculate the dimension of the global matrix
global_size = nn * 2

# Initialize the global matrix with zeros
global_matrix = np.zeros((global_size, global_size))

# Iterate over the elements and add their stiffness matrices to the global matrix
for i in range(nm):
    n1 = int(Inc[i, 0])
    n2 = int(Inc[i, 1])
    ke = lista_ke[i]

    # Determine the corresponding indices in the global matrix
    indices = [n1 * 2 - 2, n1 * 2 - 1, n2 * 2 - 2, n2 * 2 - 1]

    # Add the element stiffness matrix to the global matrix using the indices
    global_matrix[np.ix_(indices, indices)] += ke



#Aplly boundary conditions


matrix_simpl = global_matrix

matrix_simpl = np.delete(matrix_simpl, R.astype(int), axis=0)
matrix_simpl = np.delete(matrix_simpl, R.astype(int), axis=1)


forca = F
forca_simp = np.delete(F, R.astype(int), axis=0)


u2 = np.linalg.solve(matrix_simpl, forca_simp)


# print("Deslocamentos: \n",u2)
# print("Forças: \n",forca_simp)
# print("Matriz de rigidez global: \n",matrix_simpl)


u_compl = u2

for i in R.astype(int):
    u_compl = np.insert(u_compl, i, 0)

print("Deslocamentos completos: \n",u_compl)


#calcular reacoes de apoio
reacoes = np.dot(global_matrix, u_compl)
reacoes = reacoes[R.astype(int)]


# print("Reações de apoio: \n",reacoes)

#calcular tensao e deformacao
lista_deforma = []
lista_tensao = []
lista_forca_int = []
index_u = 0
for i in range(nm):
    n1 = int(Inc[i,0])
    n2 = int(Inc[i,1])        
    E = Inc[i,2]
    A = Inc[i,3]
    x1 = N[0,n1-1]
    y1 = N[1,n1-1]
    x2 = N[0,n2-1]
    y2 = N[1,n2-1]
    L = np.sqrt((x2-x1)**2 + (y2-y1)**2)
    cos = (x2-x1)/L
    sen = (y2-y1)/L
    index1 = (n1 - 1) *2
    index2 = (n2 - 1) *2
    try:
        temp_desl = np.array([u_compl[index1], u_compl[index1+1], u_compl[index2 ], u_compl[index2 +1]])
        s = E/L* np.dot(np.array([-cos, -sen, cos, sen]), temp_desl)
        e = 1/L* np.dot(np.array([-cos, -sen, cos, sen]), temp_desl)
    except:
        pass

    
    lista_deforma.append(e)
    lista_tensao.append(s)
    lista_forca_int.append(s*A)

print("Deformações: \n",lista_deforma)
print("Tensões: \n",lista_tensao)

jacob = jacobi(1000, 0.0000000000001 , matrix_simpl, forca_simp)

gauss = gaussSeidel(1000, 0.0000000000001 , matrix_simpl, forca_simp)


for i in R.astype(int):
    jacob = np.insert(jacob, i, 0)

for i in R.astype(int):
    gauss = np.insert(gauss, i, 0)



print("Deslocamentos Jacobi: \n",jacob)


print("Deslocamentos Gauss-Seidel: \n",gauss)


#compara gauss e u_compl

print("Diferença entre Gauss e u_compl: \n",gauss - u_compl)



Ft = reacoes
Ut = u_compl
Epsi = lista_deforma
Fi = lista_forca_int
Ti = lista_tensao

nome = "saidanumpy"

geraSaida(nome,Ft,Ut,Epsi,Fi,Ti)

nome = "saidajacobi"

geraSaida(nome,Ft,jacob,Epsi,Fi,Ti)

nome = "saidagauss"

geraSaida(nome,Ft,gauss,Epsi,Fi,Ti)

#reacoes 






 

