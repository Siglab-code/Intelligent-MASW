from scipy.optimize import least_squares
from scipy import optimize
import numpy as np
import disp
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from joblib import Parallel, delayed   
import os 
 



class masw():
    def __init__(self, filename, xinital,  lb, ub, Print_r=True, Savefile= True):
        self.filename = filename
        self.lb = lb
        self.ub = ub
        self.init = xinital 
        self.Print_r = Print_r
        self.Savefile = Savefile
   
        self.measure = np.loadtxt(filename) 
        self.omega = self.measure[:,0]*2*np.pi
        self.y_data =self.measure[:,1] 
        self.loss = []
        self.para = []
        self.num = (len(lb) + 1) // 4 

    
    
    def aimfuc(self, x, omega, y_data):
        vs = np.array(x[:self.num])
        mu = np.array(x[self.num:2*self.num])  
        rho = np.array(x[2 * self.num:3 * self.num]) 
        E = 2*rho*vs**2*(1+mu) 
        H1 = np.append(x[3 * self.num:], 10)     
        yt =  self.velocity(omega, E, mu, rho, H1) 
        loss = np.abs(yt - y_data)
        
        if self.Print_r: 
            plt.clf()  
            plt.scatter(omega/(2*np.pi), y_data, label= 'measurement')    
            plt.plot(omega/(2*np.pi), yt, ':', label='Prediction')    
            self.loss.append(np.sqrt(np.sum(loss**2))) 
            plt.title('TRF Algorithm, Iteration Number: %s; L2 Norm: %s' %(len(self.loss), int(self.loss[-1])))
            plt.savefig('dispersion_update.png')  
            plt.legend() 
            plt.xlabel('frequency (Hz)') 
            plt.ylabel('Phase Velocity (m/s)')
            plt.savefig('dispersion_update.png')  
            print(np.sqrt(np.sum(loss ** 2)))
        
        if self.Savefile: 
            np.savetxt('prediction.txt', yt)
            self.para.append(x)
        return loss  


    def inverse(self):
        popt = least_squares(self.aimfuc, self.init, bounds=(self.lb, self.ub), method='trf', args=(self.omega, self.y_data))
        print("inversion is done") 
        r = popt.x 
        Vs = np.array(r[:self.num])
        mu =  np.array(r[self.num:self.num*2])  
        rho =  np.array(r[self.num*2:self.num*3]) 
        H = np.array(r[self.num * 3:])
        E = 2*rho*Vs**2*2*(1+mu)
        Vp =  np.sqrt(E/((1+mu)*2*rho))  
        
        print('Shear Velocity', Vs) 
        print('Piosson Ratio', mu)  
        print('P wave velcoity:', Vp)
        print('Youngs Modulus (MPa):', E/10**6) 
        print('Density',rho)
        print('Layer Thickness', H)  

        if self.Savefile:
            np.savetxt('loss.csv', self.loss)
            np.savetxt('parameters.csv', self.para)
            
    def velocity(self, omega, E, mu, rho, H1):
        def opt(i):
            def fun(k): 
                n = len(E)
                matrix = disp.f(i, E, mu, rho, H1, k, n)
                matrix1 = np.array(matrix)
                matrix2 = np.real(matrix1)/10**11
                sign, logdet = np.linalg.slogdet(matrix2)
                return sign * np.exp(logdet)
            c_test = 10
            incre = i / c_test
            root = 0.00001
            for j in range(10**6):
                past = incre
                val1 = fun(incre)
                incre =  incre - 0.01
                val2 = fun(incre)
                if (np.real(val1) * np.real(val2) <= 0):
                    root =  optimize.brentq(fun,incre, past)      
                    break 
            return (i/root)     #give one value at a frequency
    
        def final(n):
            y = opt(omega[n])
            return y   
        
        z2 = Parallel(n_jobs=-1)(delayed(final)(i) for i in range(len(omega))) 
        return z2
 
 

 