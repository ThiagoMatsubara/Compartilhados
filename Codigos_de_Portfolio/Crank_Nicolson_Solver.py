#Algoritmo de Thomas para solucionar o método de Crank - Nicolson
#informações
r'''
Equação implícita 4.40 - John Anderson

\frac{(T_{i}^{n+1}-T_{i}^{n})}{\Delta t} = \alpha \cdot \frac{1}{2} \cdot \frac{[(T_{i+1}^{n+1}-T_{i+1}^{n})+(2T_{i}^{n+1}-2T_{i}^{n})+(T_{i-1}^{n+1}-T_{i-1}^{n})]}{(\Delta x)^2}

Separando do lado esquerdo os termos não disponíveis: 

\frac{\alpha \Delta t}{2 (\Delta x)^2}T^{n+1}_{i+1} - [1+\frac{\alpha \Delta t}{2 (\Delta x)^2}]T^{n+1}_{i} + \frac{\alpha \Delta t}{2 (\Delta x)^2}T^{n+1}_{i-1} = -T^{n}_{i} - \frac{\alpha \Delta t}{2 (\Delta x)^2}(T^{n}_{i+1}-2T^{n}_{i}+T^{n}_{i-1})

simplificando a nomenclatura:


A = \frac{\alpha \Delta t}{2 (\Delta x)^2}

\\ \\

B = [1+\frac{\alpha \Delta t}{2 (\Delta x)^2}]

\\ \\

K_{i} = -T^{n}_{i} - \frac{\alpha \Delta t}{2 (\Delta x)^2}(T^{n}_{i+1}-2T^{n}_{i}+T^{n}_{i-1})


Aplicando no grid de pontos de 2 a 6, ou seja nos valores de i, como se sabe têm-se os valores de T1 e T7 (no caso do grid de 7 pontos, 1 e 7 são condições de contorno)

chega-se a matriz: 

[ [-B, A, 0, 0, 0],
  [A, -B, A, 0, 0],
  [0, A, -B, A, 0],
  [0, 0, A, -B, A],
  [0, 0, 0, A, -B]
]

*

[T2,T3,T4,T5,T6]

=

[K2',K3,K4,K5,K6']

onde K2' = K2 -AT1
onde K6' = K6 -AT7

Essa matriz trigonal é resolvida pelo algoritmo de Thomas

Informações: 
Diferenciação usada principalmente para resolver problemas governados por equações parabólicas

Algoritmo de solução:

Não mechemos na primeira nem ultima linha
Da segunda até a penúltima linha: 
	Eliminados o primeiro termo da esquerda para direita diferente de zero
	Trocamos o termo da diagonal da matriz e o termo do outro lado da igualdade

Exemplo:
	
[
[B1,C1,0,0,0,0,0],
[A2,B2,C2,0,0,0,0],
[0,A3,B3,C3,0,0,0],
[0,0,A4,B4,C4,0,0],
[0,0,0,A5,B5,C5,0],
[0,0,0,0,A6,B6,C6],
[0,0,0,0,0, 0, B7],
]

*

[T1,T2,T3,T4,T5,T6,T7]

=

[K1,K2,K3,K4,K5,K6, K7]

PÓS OPERAÇÃO:

[
[B1,C1,0,0,0,0,0],
[0,B2',C2,0,0,0,0],
[0,0,B3',C3,0,0,0],
[0,0,0,B4',CA4,0,0],
[0,0,0,0,B5',C5,0],
[0,0,0,0,0,B6',C6],
[0,0,0,0,0, 0, B7'],
]

*

[T1,T2,T3,T4,T5,T6,T7]

=


[K1,K2',K3',K4',K5',K6', K7]

Sendo que:

B2' = B2 - C1*A2/B1
B3' = B3'- C2*A3/B2'
...
Bi' = Bi' - C_(i-1)*Ai/B_(i-1)'


e 
[K1,K2',K3',K4',K5',K6', K7] ->

K2' = K2 - K1*A2/B1
K3' = K3 - K2*A3/B2'
...
Ki' = Ki - K_(i-1)*Ai/B'_(i-1)

Substituindo todos esses termos, é possível encontrar o ultimo termo do vetor T (das temperaturas neste caso)
Tn = Kn'/Bn' e, consequentemente todos os outros podem ser encontrados de baixo para cima

T_(n-1) = [K_(n-1)' - C_(n-1)*T_(n)]/B_(n-1)' 

ou

T_(n) = [K_(n)' - C_(n)*T_(n+1)]/B_(n)'


** Note que a nomenclatura mudou, na matriz mostrada antes do tópico algoritmo de solução
possuía duas linhas a menos pois K7 e K1 são conhecidos, aqui ele é "escrita por extenso"

DESVANTAGENS: 
A desvantagem desse método é que as condições de contorno precisam estar no primeiro
e ultimo nó no caso 1D, e em todas as paredes no caso 2D
Além disso, o caso 2D é complexo de resolver matemáticamente.
'''

import numpy as np
import matplotlib.pyplot as plt
#Entrada de parâmetros:
#Quanto maior o numero de nós, menor deve ser o passo de tempo para manter a estabilidade numérica
alpha = 80
L = 10
nodes = 1000
delta_x = L/(nodes-1)
delta_t = 0.001
t = 3
T_i = 15

A = alpha*delta_t/(2*(delta_x**2))
B = np.zeros(nodes) + (1+2*A)*(-1)
K = np.zeros(nodes)
T = np.zeros(nodes) + T_i

#Boundary conditions, tem que ser em 0 e -1!
T[0] = 65
T[-1] = 100


#


# Preciso solucionar a equação:
# A*T_(i+1)^(n+1) - B*T_(i)^(n+1)+A*T_(i-1)^(i+1) = Ki
# Somente a temepratura esta no dominio do tempo e do espaço,
# K e B estão apenas no dominio do espaço
# A matrix que estamos solucionando é a da equação 4.51, claro, com o 
# numero de nós adaptados para o caso.

fig = plt.figure()#figsize = (12,12))

axis1 = plt.subplot(2, 2, 1)
pcm = axis1.pcolormesh([T], cmap=plt.cm.jet, vmin=min(T), vmax=max(T))
axis1.set_ylim([-5,5])
plt.colorbar(pcm,ax=axis1)




print(f'B antes de qualquer modificação: {B}')
time =0
while time < t:
    T_ant = np.copy(T)
 
    for i in range(1,nodes-1):
        K[i] = -T_ant[i] - A*(T_ant[i+1]-2*T_ant[i]+T_ant[i-1])

    #das condições de contorno
    #Equivale ao K2' e o K6' do livro
    K[1] = K[1] - A*T[0]
    K[nodes-2]= K[nodes-2] - A*T[-1]

    #Agora, Resolvendo a matriz:
    K_ant = np.copy(K)
    B_ant = np.copy(B)
    for i in range(2, nodes-1):
        
        B_ant[i] = B_ant[i] - (A**2)/B_ant[i-1]
        K_ant[i] = K_ant[i] - (K_ant[i-1]*A)/B_ant[i-1]

    #A ultima linha da matriz é a primeira a ser resolvida
    T[nodes-2] = K_ant[nodes-2]/B_ant[nodes-2]

    #De baixo para cima: 
    for i in range(nodes-3,0,-1):
        T[i] = (K_ant[i] - A*T[i+1])/B_ant[i]

    time += delta_t

    pcm.set_array([T])
    axis1.set_title("Distribution at t: {:.3f}s".format(time))

    plt.pause(0.01)
plt.show()


print(f"Perfil de Temperatura em t={time}s:")
print(T)
plt.plot(np.linspace(0, L, nodes), T, marker='o', label=f't = {time}s')
plt.xlabel('Posição (mm)')
plt.ylabel('Temperatura (°C)')
plt.title('Crank-Nicolson (Algoritmo de Thomas)')
plt.grid(True)
plt.show()