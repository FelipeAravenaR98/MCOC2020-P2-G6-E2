import numpy as np

g = 9.81 #kg*m/s^2


class Barra(object):

    """Constructor para una barra"""
    def __init__(self, ni, nj, R, t, E, ρ, σy):
        super(Barra, self).__init__()
        self.ni = ni
        self.nj = nj
        self.R = R
        self.t = t
        self.E = E
        self.ρ = ρ
        self.σy = σy

    def obtener_conectividad(self):
        return [self.ni, self.nj]

    def calcular_area(self):
        A = np.pi*(self.R**2) - np.pi*((self.R-self.t)**2)
        return A

    def calcular_largo(self, reticulado):
        """Devuelve el largo de la barra. 
        ret: instancia de objeto tipo reticulado
        """
        xi = reticulado.obtener_coordenada_nodal(self.ni)
        xj = reticulado.obtener_coordenada_nodal(self.nj)
        dij = xi-xj
        return np.sqrt(np.dot(dij,dij))

    def calcular_peso(self, reticulado):
        """Devuelve el largo de la barra. 
        ret: instancia de objeto tipo reticulado
        """
        L = self.calcular_largo(reticulado)
        A = self.calcular_area()
        return self.ρ * A * L * g






    def obtener_rigidez(self, ret):
        """Devuelve la rigidez ke del elemento. Arreglo numpy de (4x4)
        ret: instancia de objeto tipo reticulado
        """
        L = self.calcular_largo(ret)
        A = self.calcular_area()
        
        k = self.E * A / L 
        
        
        xi=ret.xyz[self.ni,0]
        yi=ret.xyz[self.ni,1]
        
        xj=ret.xyz[self.nj,0]
        yj=ret.xyz[self.nj,1]
        
        #Segun el video T teta es delta definido asi:
        menos_cos = (xi-xj) / L
        menos_sen = (yi-yj) / L
        cos       = (xj-xi) / L
        sen       = (yj-yi) / L

        T0= np.array([[menos_cos], [menos_sen], [cos], [sen]])   
        
        ke = T0 @ T0.T * k
        #implementar

        return ke

    def obtener_vector_de_cargas(self, ret):
        """Devuelve el vector de cargas nodales fe del elemento. Vector numpy de (4x1)
        ret: instancia de objeto tipo reticulado
        """
        W = self.calcular_peso(ret)
        
        fe = np.array([[0],[-1],[0],[-1]]) * (W / 2)
        #Implementar

        return fe


    def obtener_fuerza(self, ret):
        """Devuelve la fuerza se que debe resistir la barra. Un escalar tipo double. 
        ret: instancia de objeto tipo reticulado
        """
        #Implementar
        A=self.calcular_area()
        L=self.calcular_largo(ret)        
        ni=self.ni
        nj=self.nj
        ue=np.array([ret.u[ni*2], ret.u[((2*ni)+1)], ret.u[2*nj], ret.u[((2*nj)+1)]])       
        xi=ret.xyz[ni,0]
        yi=ret.xyz[ni,1]
        
        xj=ret.xyz[nj,0]
        yj=ret.xyz[nj,1]
        
        menos_cos = (xi-xj) / L
        menos_sen = (yi-yj) / L
        cos       = (xj-xi) / L
        sen       = (yj-yi) / L

        T0= np.array([[menos_cos], [menos_sen], [cos], [sen]])
        delta = T0.T @ ue
        
        se = ((A*self.E) / L )* delta

        return se























