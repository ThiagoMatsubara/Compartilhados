import numpy as np
import matplotlib.pyplot as plt

class Conduction:
    """
    MANUAL DO CASO
    A equação de condução de calor unidimensional é dada por:
    Del(T)/del(t) = alpha* del²(T)/del(x)²


    Aqui é importante saber que:
    Serie de Taylor
    f(x+delta(x)) = f(x) + del(f)/del(x) * delta(x)/1! + del²(f)/del(x)² * delta(x)²/2! + del³(f)/del(x)³ + delta(x)³/3! + ...
    Sendo que a é o proprio valor de x, (x_Final - a=x)=delta(x)

    A diferencial parcial pode ser aproximada pelo método das diferenças finitas, de três maneiras:
    Foward difference
    Rearward difference
    Central difference

    time
    ^
    | n + 1     .       .       .
    |  n        .       .       .
    | n - 1     .       .       .
    |           i-1     i       i+1
    |           ----------------------> x

    Usando a foward difference para a hiperbólica: 
    (delt(T)/del(t))_n,i = (T_n+1,i - T_n,i)/delta(t) - (del²(T)/delt(t))_n,i * delta(t)/2 +...
    i is the running index in the x direction
    n is the running index in the t direction (avanço no tempo)

    Usando as diferenças centrais para a parabolica:
    del²(T)/del(x)² = (T_i+1,j - 2T_i,j + T_i-1,j)/(delta(x))² + O(delta(X))²


    Os Gráficos auxiliares representam:
    Residuos lineares: Máxima diferença entre a solução atual e anterior para cada variável (se descendo, convergindo)
    Gráfico da temperatura média da barra ao longo do tempo
    Gráfico da temepratura ao longo da barra em cada iteração

    MANUAL DE USO:

        from conduction import Conduction


        dados_1D = Conduction(alpha=80,L=40,t=2.5,nodes=40, dim = True)
        cond_1D = {(0,10):80, (-1,None):80}
        dados_1D.def_boundary_conditions(cond_1D)
        dados_1D.conducao_1D()


        dados_2D = Conduction(alpha=80,L=40,t=2.5,nodes=40, dim = False)
        cond_2D = {(-1, None, 0, None):150, (17,23,17,23):150, (0, 1, 0, None):150, (0,None,0,1):150, (0,None,-1,None):150}
        dados_2D.def_boundary_conditions(cond_2D)
        dados_2D.conducao_2D()
    """

    def __init__(self,alpha,L=10.0,t=5.0,nodes=10,t_i=20.0,dim = True):

        self.alpha = alpha
        self.L = L
        self.t = t
        self.nodes = nodes
        self.t_i = t_i
        self.dim = dim


        self.dx = L/(nodes-1)
        self.dy = L/(nodes-1)

        if dim == True:
            self.l = np.zeros(nodes) + t_i
            self.dt = 0.5 * self.dx**2 / alpha
        elif dim == False:
            self.l = np.zeros((nodes, nodes)) + t_i
            self.dt= min(self.dx**2/(4*self.alpha),self.dy**2/(4*self.alpha))



    def __str__(self):
        return f'alpha: {self.alpha:.2f}\nL: {self.L:.2f}\nt: {self.t:.2f}\nnodes: {self.nodes:.2f}'

    def def_boundary_conditions(self, cond = None):

        '''
        Você pode não passar nenhuma condição, ou passar alguma. Por exemplo:
        minhas_conds = {(0, 3): 80, (37, 40): 100} para o caso 1D
        minhas_conds = {(-1,None,1,None):80, (3,6,3,6):150}

        x = Conduction(alpha=250, L=500, t=15, nodes=40)
        x.def_boundary_conditions(condicoes=minhas_conds)


        '''

        if cond:

            if self.dim == True:
                for intervalo, temperatura in cond.items():
                    a,b = intervalo
                    self.l[a:b]=temperatura
                print('condições aplicadas!')

            if self.dim == False:
                for intervalo, temperatura in cond.items():
                    a,b,c,d = intervalo
                    self.l[a:b,c:d] = temperatura
                print('condições aplicadas!')

        else:
            print('Defina as condições de contorno:\n')
            resposta = 'S'
            while resposta == 'S':
                
                if self.dim == True:
                    print('Digite um valor para o intervalo da condição de contorno na barra do tipo [a:b]=c')
                    print('Exemplo: [0:2]=80. Significa que a temperatura da barra sera 80 nos primeiros 3 nós da barra')
                    a = input('a: ').strip()
                    b = input('b: ').strip()
                    c = float(input('c: '))

                    a = int(a) if a else None
                    b = int(b) if b else None
                    if a is not None and b is None and b =='':
                        self.l[a:] = c
                    
                    if b is not None and a is None and a =='':
                        self.l[:b] = c
                    
                    else:
                        self.l[a:b]=c
                
                elif self.dim == False:
                    print('Digite um valor para o intervalo da condição de contorno na barra do tipo [a:b, c:d]=e')
                    print('Exemplo: [0:2, 0:10]=80. Significa que a temperatura da placa sera 80 nos tres primeiros nós verticais, ' \
                    'e nos 11 primeiros nós horizontais')

                    a = input('a: ').strip()
                    b = input('b: ').strip()
                    c = input('c: ').strip()
                    d = input('d: ').strip()
                    e = float(input('e: '))

                    a = int(a) if a else 0
                    b = int(b) if b else None
                    c = int(a) if c else 0
                    d = int(b) if d else None
                    
                else:
                    self.l[a:b, c:d]=e
                
                
                v=0
                while v==0:
                    print('Deseja inserir uma condição de contorno?')
                    resposta = str(input('S/N: ')) 
                    if resposta == 'S':
                        v=1
                    if resposta == "N":
                        v=1
                        break
                    if resposta != "N" and resposta != "S":
                        print("digite um comando válido!")
            
    def conducao_1D(self):

        #Simulação:
        fig = plt.figure()#figsize = (12,12))

        axis1 = plt.subplot(2, 2, 1)
        pcm = axis1.pcolormesh([self.l], cmap=plt.cm.jet, vmin=0, vmax=100)
        axis1.set_ylim([-5,5])
        plt.colorbar(pcm,ax=axis1)


        axis2 = plt.subplot(2, 2, 2)
        axis2.set_xlabel('t (s)')
        axis2.grid()

        axis3 = plt.subplot(2, 2, 3)
        axis3.set_xlabel('t (s)')
        axis3.grid()

        axis4 = plt.subplot(2,2,4)
        axis4.set_title('Temeperaturas em cada iteração')
        axis4.set_xlabel('L(mm)')
        axis4.set_ylabel('T(°C)')
        axis4.grid()

        plt.tight_layout()

        #fig, ((axis1,axis2), axis3) = plt.subplots(2,1, figsize = (10,10))

        contador = 0
        time_graf2 = []
        avarege_l=[]

        residuals=np.zeros(self.nodes)
        log_residuo = []
        log_temps = []

        while contador < self.t: 
            w = self.l.copy()
            r = residuals.copy()

            time_graf2.append(contador)
            avarege_l.append(np.average(self.l))
            log_residuo.append(max(residuals))
            log_temps.append(self.l)

            for i in range(1, self.nodes -1):
                #T_(n+1),i = T_n,i + alpha * delta(t)/(delta(x)²) * (T_t,i+1 -2T_t,i + T_t,i-1)
                self.l[i] = w[i] + self.alpha * self.dt * (w[i+1] - 2 * w[i] + w[i -1])/(self.dx**2)
                residuals[i] = self.l[i] - w[i]

            contador += self.dt


            #print(f"time: {contador:.3f}s, Temperatura: {np.average(l):.2f} ")
            #print(f'Temperatura na parede: {l[0]}')

            pcm.set_array([self.l])
            axis1.set_title("Distribution at t: {:.3f}s".format(contador))

            axis2.plot(time_graf2,avarege_l, color = 'black', linestyle = ':')
            axis2.set_title("Average Temperature Evolution: {:.2f}".format(np.average(self.l)))

            axis3.plot(time_graf2,log_residuo, color = 'red', linestyle='--')
            axis3.set_title(f'Residuals: {max(residuals):.2f}')

            axis4.plot(np.linspace(1,self.nodes,self.nodes),log_temps[0])
            
            
            
            plt.pause(0.01)

        plt.show()

    def conducao_2D(self):

        '''
        # A equação de condução de calor bidimensional é dada por:
        # Del(T)/del(t) = alpha* (del²(T)/del(x)² + del²(T)/del(y)²)


        #Aqui é importante saber que:

        #Serie de Taylor:
        #f(x+delta(x)) = f(x) + del(f)/del(x) * delta(x)/1! + del²(f)/del(x)² * delta(x)²/2! + del³(f)/del(x)³ + delta(x)³/3! + ...

        #Sendo que a é o proprio valor de x, (x_Final - a=x)=delta(x)

        #A diferencial parcial pode ser aproximada pelo método das diferenças finitas, de três maneiras:
        # Foward difference
        # Rearward difference
        # Central difference


        #Usando a foward difference para a hiperbólica: 
        #(delt(T)/del(t))_n,i = (T_n+1,i,j - T_n,i,j)/delta(t) - (del²(T)/delt(t))_n,i * delta(t)/2 +...

        # i is the running index in the x direction
        # j is the running index in the y direction
        # n is the running index in the t direction (avanço no tempo)

        #Usando as diferenças centrais para a parabolica:
        #del²(T)/del(x)² = (T_i+1,j - 2T_i,j + T_i-1,j)/(delta(x))² + O(delta(X))²


        # A lembrar, ou talvez, saber pela primeira vez
        # Para declarar as condições de contorno, pode-se navegar pelas paredes da seguinte maneira; 

        #           [-1,:]
        #         _ _ _ _ _ _
        #       |            |
        # [:,1] |            |  [:,-1]
        #       | _ _ _ _ _ _|
        #           [0,:]
        
        ou seja, linhas depois colunas
        '''

        #Simulação:
        fig = plt.figure()#figsize = (12,12))

        axis1 = plt.subplot(2, 2, 1)
        pcm = axis1.pcolormesh(self.l, cmap=plt.cm.jet, vmin=np.min(self.l), vmax=np.max(self.l))
        plt.colorbar(pcm,ax=axis1)


        axis2 = plt.subplot(2, 2, 2)
        axis2.set_xlabel('t (s)')
        axis2.grid()

        axis3 = plt.subplot(2, 1, 2)
        axis3.set_xlabel('t (s)')
        axis3.grid()

        plt.tight_layout()

        #fig, ((axis1,axis2), axis3) = plt.subplots(2,1, figsize = (10,10))

        contador = 0
        time_graf2 = []
        avarege_l=[]

        residuals=np.zeros((self.nodes, self.nodes))
        log_residuo = []
        log_temps = []

        while contador < self.t: 
            w = self.l.copy()
            r = residuals.copy()

            time_graf2.append(contador)
            avarege_l.append(np.average(self.l))
            log_residuo.append(np.max(np.abs(residuals)))
            log_temps.append(self.l)

            for i in range(1, self.nodes -1):
                for j in range(1,self.nodes -1):
                    self.l[i,j] = w[i,j] + self.alpha * self.dt * ( (w[i+1,j] - 2 * w[i,j] + w[i -1,j])/(self.dx**2) +
                                            (w[i,j+1] - 2 * w[i,j] + w[i,j-1])/(self.dy**2))
                
                    residuals[i,j] = self.l[i,j] - w[i,j]

            contador += self.dt


            #print(f"time: {contador:.3f}s, Temperatura: {np.average(l):.2f} ")
            #print(f'Temperatura na parede: {l[0]}')

            pcm.set_array(self.l)
            axis1.set_title("Distribution at t: {:.3f}s".format(contador))

            axis2.plot(time_graf2,avarege_l, color = 'black', linestyle = ':')
            axis2.set_title("Average Temperature Evolution: {:.2f}".format(np.average(self.l)))

            axis3.plot(time_graf2,log_residuo, color = 'red', linestyle='--')
            axis3.set_title(f'Residuals: {np.max(np.abs(residuals)):.2f}')
            
            plt.pause(0.01)

        plt.show()

