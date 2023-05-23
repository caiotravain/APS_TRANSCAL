from funcoesTermosol import *
import numpy as np
import matplotlib.pyplot as plt

def matprint(mat, fmt="g"):
    col_maxes = [max([len(("{:"+fmt+"}").format(x)) for x in col]) for col in mat.T]
    for x in mat:
        for i, y in enumerate(x):
            print(("{:"+str(col_maxes[i])+fmt+"}").format(y), end="  ")
        print("")
  

[nn,N,nm,Inc,nc,F,nr,R] = importa('entrada.xlsx')

#plota(N,Inc)

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


global_matrix_2 = global_matrix

global_matrix_2 = np.delete(global_matrix_2, R.astype(int), axis=0)
global_matrix_2 = np.delete(global_matrix_2, R.astype(int), axis=1)


forca = F
forca2 = np.delete(F, R.astype(int), axis=0)


u2 = np.linalg.solve(global_matrix_2, forca2)



u = np.zeros((nn*2,1))

for i in range(len(u)):
    if i not in R:
        u[i] += u2[i%len(u2)]




print("Deslocamentos: \n",u)





#reacoes 






    
#print(lista_ke[0])
# lista_kg[0:4,0:4] += lista_ke[0]
# lista_kg[0:2,0:2] += lista_ke[2][2:4,2:4]
# lista_kg[2:4,2:4] += lista_ke[2][0:2,0:2]



