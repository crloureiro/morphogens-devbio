import numpy as np
            
class Diffusion1D:
    
    def __init__(self,N):
        self.N = N
        self.a = np.zeros(N)
        self.b = np.zeros(N)
        self.c = np.zeros(N)
        self.d = np.zeros(N)
        self.e = np.zeros(N)
        self.f = np.zeros(N)
        self.s = np.zeros(N)
        self.alpha = np.zeros(N)
        self.beta = np.zeros(N)
        self.gamma = np.zeros(N)
        self.x = np.zeros(N)
    
    def config(self,dt,coef,type0,lim0,type1,lim1,S):
        N = self.N
        dx = 1.0/(N-1)
        self.dt = dt
        v = 2*dt/(dx*dx)
        i1 = 0
        for k in range(len(coef)): # coefficients de diffusion
            x = coef[k][0]
            diff = coef[k][1]
            i2 = int(x*N)
            for i in range(i1,i2):
                self.alpha[i] = diff*v
            i1 = i2
        for i in range(1,N-1): # Cranck-Nicolson
            self.a[i] = -self.alpha[i]/2
            self.b[i] = (1+self.alpha[i])
            self.c[i] = self.a[i]
            self.d[i] = self.alpha[i]/2
            self.e[i] = 1-self.alpha[i]
            self.f[i] = self.alpha[i]/2
            self.s[i] = dt*S[i]
        if len(coef)>1: # frontieres entre deux milieux differents
            for k in range(len(coef)-1):
                i = int(coef[k][0]*N)
                self.a[i] = -coef[k][1]
                self.b[i] = coef[k][1]+coef[k+1][1]
                self.c[i] = -coef[k+1][1]
                self.d[i] = 0
                self.e[i] = 0
                self.f[i] = 0
                self.s[i] = 0
        if type0=="dirichlet":
            self.b[0] = 1.0
            self.c[0] = 0.0
            self.e[0] = 0
            self.f[0] = 0
            self.s[0] = lim0
        if type1=="dirichlet":
            self.a[N-1] = 0
            self.b[N-1] = 1.0
            self.d[N-1] = 0
            self.e[N-1] = 0
            self.s[N-1] = lim1
        if type0=="neumann":
            self.b[0] = 1.0+self.alpha[0]
            self.c[0] = -self.alpha[0]
            self.e[0] = 1.0-self.alpha[0]
            self.f[0] = self.alpha[0]
            self.s[0] = -2*self.alpha[0]*lim0*dx+dt*S[0]
        if type1=="neumann":
            self.a[N-1] = -self.alpha[N-1]
            self.b[N-1] = 1.0+self.alpha[N-1]
            self.d[N-1] = self.alpha[N-1]
            self.e[N-1] = 1.0-self.alpha[N-1]
            self.s[N-1] = 2*self.alpha[N-1]*dx*lim1+dt*S[N-1]
        
        self.beta[0] = self.b[0]
        self.gamma[0] = self.c[0]/self.beta[0]
        for k in range(1,N):
            self.beta[k] = self.b[k]-self.a[k]*self.gamma[k-1]
            if self.beta[k]==0:
                raise Exception("Impossible de resoudre le systeme ligne "+str(k))
            self.gamma[k] = self.c[k]/self.beta[k]
           
    def iterations(self,U0,ti,tf):
        if U0.size!=self.N:
            raise Exception("Taille de U incorrecte")
        U = U0.copy()
        t = ti
        while t<tf:
            self.x[0] = (self.e[0]*U[0]+self.f[0]*U[1]+self.s[0])/self.beta[0]
            for k in range(1,self.N-1):
                self.x[k] = (self.d[k]*U[k-1]+self.e[k]*U[k]+self.f[k]*U[k+1]+self.s[k]-self.a[k]*self.x[k-1])/self.beta[k]
            k = self.N-1
            self.x[k] = (self.d[k]*U[k-1]+self.e[k]*U[k]+self.s[k]-self.a[k]*self.x[k-1])/self.beta[k]
            U[self.N-1] = self.x[self.N-1]
            for k in range(self.N-2,-1,-1):
                U[k] = self.x[k]-self.gamma[k]*U[k+1]
            t += self.dt
        return [U,t]
            