import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt

c=299792.458 #speed of light in km/s 
H0=70.0     #Hubble Constant in km/s.Mpc

omega_lam=0.70
omega_m=1.0-omega_lam

def DH(H0):    #Hubble Distance  
    
    return c/H0
    
def Integrand(z,omega_m,omega_lam):  #The E(z) function given in D.W. Hogg (2000) arXiv:astro-ph/9905116
    
    #return (omega_m*(1+z)**3+omega_lam)**(-0.5)
    return ((1+z)**2*(1+omega_m*z)-z*(2+z)*omega_lam)**(-0.5)



def DC(z):  #The comoving distnace
    
    #comovdist=DH(H0)*quad(Integrand, 0, z, args=(omega_m,omega_lam))[0] 
    comovdist=quad(Integrand, 0, z, args=(omega_m,omega_lam))[0]
    return comovdist

def DL(z): #The luminosity distance
    
    return (1+z)*DC(z)
    
def mu(z):  #The distance modulus
    
    return 5*np.log10(DL(z))+25
    


z=0.36#00342119890000

ans=5*np.log10(DL(z))

            
#DC=np.vectorize(DC) 
#z = np.linspace(0,1.5,100) 
#dcomov = DC(z)

#DL=np.vectorize(DL) 
#z = np.linspace(0,5.0,100) 
#lumdist = DL(z)

#mu=np.vectorize(mu) 
#z = np.linspace(0.01,1.5,100) 
#distmodul = mu(z)


a,b,c=np.genfromtxt("Union2.1_z_dm_err.txt",unpack=True)
#plt.plot(a,b)
#plt.errorbar(a,b,c, fmt= 'bo')

#plt.plot(z,distmodul, linewidth=3, color='black')

#plt.savefig('/Users/behnam/work/python_SN/mu_vs_z.eps')




plt.show()
