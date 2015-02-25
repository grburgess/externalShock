from matplotlib.pyplot import *
from numpy import *
from scipy.integrate import quad, quadrature


class extShock_adi(object):
    
    def __init__(self,z=1.,E0=1.E54,Gamma0=300.,eta=0.,j=2.,n0=1.E2,q=1.E-3):
        
        
        from astropy import cosmology
        from astropy.constants import c


        cosmo=cosmology.default_cosmology.get()
        
        self.c      = c.cgs.value
        
        self.mp     = 0.00150327741 # erg   #1.67E-26 #keV?? WHAT THE FUCK ARE THESE UNITS????


        self.j = j

        g = (j+1.-eta)*0.5
        
        self.eta    = eta
        self.Gamma0 = Gamma0
        self.g      = g
        self.n0     = n0
        self.z      = z
        self.dL = (cosmo.luminosity_distance(self.z)).cgs.value
        self.E0     = E0
        
        self.nu     = 4./3.
        self.delta  = .2
        self.q      = q
        
        
        self.xd = 1.8E16*((1.+self.j-self.eta)*(self.E0/1.E54)/((self.n0/100.)*(self.Gamma0/300.)))**(1./3.)
        
        self.td = 6.7*(1.+self.z)*((1.+self.j-self.eta)*(self.E0/(1.E54))/((self.n0/100.)*(self.Gamma0/300.)**8.))**(1./3.)

        A0 = 4*pi*self.xd**(2.)
        self.A0 = A0

        #Here there is a modification accounting for the area scaling: j
        self.Pi0  = ((2.*self.g - self.j+1. +self.eta)*self.mp*self.c*self.Gamma0**4.*self.n0*A0)
        self.Pi0 /= (2.*self.g*(1./self.nu +1./self.delta)*(1.+self.z)**2.)
#        print self.Pi0
        
        self.eE0 = 3.E-8 * self.n0**(.5)*self.q*self.Gamma0**4. /(1.+self.z)
        self.sedovRadius = self.Gamma0**(2./3.)*self.xd
    
    
    def Gamma(self,x):
    
        if x <1.:
            return self.Gamma0
    
        elif (1<=x) and (x<=self.Gamma0**(1./self.g)):
        
            return self.Gamma0*x**(-self.g)




    def ep(self, t):
    
        val = self.eE0*(self.Gamma(self.X(t))/self.Gamma0)**4. * (self.X(t))**(-self.eta/2.)
        #x=self.X(t)
        #if x>=0 and x<1:
        #    val  = x**(-self.eta/2.)
        #elif x>=1 and x<self.Gamma0**(1./self.g):
        #    val = x**(-self.eta/2.-4.*self.g)
       
        #else:
        #    val = 0.
            
        #val*=self.eE0
        
        return val





    def P_p(self, t):
    
    
        #val = self.Pi0*(self.Gamma(self.X(t))/self.Gamma0)**4. * (self.X(t)/self.xd)**(2.-self.eta)
        x=self.X(t)
        if x>=0 and x<1:
            val  = x**(self.j-self.eta)
        elif x>=1 and x<self.Gamma0**(1./self.g):
            val = x**(self.j-self.eta-4.*self.g)
            
        else:
            val = 0
        val*=self.Pi0
    
    
        return val


    def P(self,ene,t):
    
        top = (1.+self.nu/self.delta)*self.P_p(t)
        bottom = (ene/self.ep(t))**(-self.nu) + (self.nu/self.delta)*(ene/self.ep(t))**(self.delta)
    
        val = top/bottom
        val/=(4*pi*self.dL**2)

        return val ## BLAHHHHH
    
    
    def X(self,t):
    
    
        frac = t/self.td 
        
        test = (self.td/(2. * self.g + 1.)) * (self.Gamma0**(2.+1./self.g) + 2.*self.g)
        if t<self.td:
        
            val = frac

            if val >=0. and val <1:
                return val
        
        
        elif self.td<=t and t<=test:
        
            val = ((2.*self.g+1.)*frac - 2.*self.g)**(1./(2.*self.g+1.))
            
            if val>=1 and val<= self.Gamma0**(1./self.g):
                
                return val
        else:
            
            return 0.
        
        
    
    def PlotGamma(self,tdStart=.1,tdStop=1E1,nT=10):
    
        tGrid = logspace(log10(tdStart*self.td),log10(tdStop*self.td),nT)
        
        xs = array(map(lambda t: self.X(t),tGrid))
        
        gammas = array(map(lambda x: self.Gamma(x), xs ))
        
        
        loglog(tGrid,gammas)
        
       
        
    def GetRadii(self,tdStart=.1,tdStop=1E1,nT=10): 
        tGrid = logspace(log10(tdStart*self.td),log10(tdStop*self.td),nT)
        vals = array(map(lambda t: self.X(t),tGrid))*self.xd
        
        return vals
    
    
    def GetGamma(self,tdStart=.1,tdStop=1E1,nT=10):
        
        tGrid = logspace(log10(tdStart*self.td),log10(tdStop*self.td),nT)
        
        xs = array(map(lambda t: self.X(t),tGrid))
        
        gammas = array(map(lambda x: self.Gamma(x), xs ))
        
        return gammas
        
        
    def PlotSpecEvo(self,tdStart=.1,tdStop=1E1,nT=10):
        
        tGrid = logspace(log10(tdStart*self.td),log10(tdStop*self.td),nT)
        eGrid = logspace(-6.,6,100.)
        for t in tGrid:
            spec = array(map(lambda e: es.P(e,t), eGrid ))
            loglog(eGrid,eGrid**-0.*spec)
        
        
    def PlotEpEvo(self,tdStart=.1,tdStop=1E1,nT=10,color='red'):
        
        
        tGrid = logspace(log10(tdStart*self.td),log10(tdStop*self.td),nT)
            
        eps = array(map(lambda t: 511.*self.ep(t)  ,tGrid))
        
        
        loglog(tGrid,eps,color='color')
      
    def PlotPeakFluxEp(self,tdStart=.1,tdStop=1E1,nT=10,color='red'):
        
        
        tGrid = logspace(log10(tdStart*self.td),log10(tdStop*self.td),nT)
        eps = array(map(lambda t: 511.*self.ep(t)  ,tGrid))
        ps = array(map(lambda t: self.P_p(t)  ,tGrid))
        
        loglog(eps,ps,color=color)
         
    def PlotLightCurve_vFv(self,tdStart=.1,tdStop=1E1,nT=100,ene=1. ):
        
        tGrid = logspace(log10(tdStart*self.td),log10(tdStop*self.td),nT)
        
        lc = array(map(lambda t: self.P(ene,t) , tGrid))
        
        plot(tGrid,lc)

    def PlotLightCurve(self,tdStart=.1,tdStop=1E1,nT=100,ene=1. ):
        
        tGrid = logspace(log10(tdStart*self.td),log10(tdStop*self.td),nT)
        
        lc = array(map(lambda t: self.P(ene,t)/(511.*ene) , tGrid))
        
        plot(tGrid,lc)

            
    def IntegralFlux(self,t,emin=10.,emax=40000.):


        emin = emin/511.
        emax = emax/511.


        
