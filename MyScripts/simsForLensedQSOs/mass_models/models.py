import numpy as np


class _MassModel:
    """
    A generic mass model class, including a method to calculate deflection
        angles for power law models. This probably needs to be re-worked.
    """
    def align_coords(self,xin,yin,revert=False):
        theta = self.theta-np.pi/2.
        ctheta = np.cos(theta)
        stheta = np.sin(theta)
        if revert:
            X = xin*ctheta-yin*stheta
            Y = yin*ctheta+xin*stheta
            x = X+self.x
            y = Y+self.y
            return x,y
        X = xin-self.x
        Y = yin-self.y
        x = X*ctheta+Y*stheta
        y = Y*ctheta-X*stheta
        return x,y


    def get_deflection_angles(self,x,y,**args):
        """
        Powerlaw deflection angles following Barkana
        """
        x = x.ravel()
        y = y.ravel()
        q = self.q
        x,y = self.align_coords(x,y)
        x2 = x**2
        y2 = y**2
        rho0 = (x2+y2/q**2)**0.5
        ell = (1.-q**2)**0.5

        rho = np.linspace(1e-11,1.,self.nsteps)
        rho = rho.repeat(x.size).reshape((self.nsteps,x.size))
        rho *= rho0

        rhoell = (rho*ell)**2
        d = ((rhoell+y2-x2)**2 + 4*x2*y2)**0.5
        w = ((d+x2+y2+rhoell)/(d+x2+y2-rhoell))**0.5

        kappa = self.kappa(rho,**args)
        integrand = rho*kappa*w/(x2+y2*w**4)
        xmap = 2.*x*q*integrand.sum(0)*rho0/self.nsteps
        ymap = 2.*y*q*(integrand*w**2).sum(0)*rho0/self.nsteps
        xmap[rho0==0] = 0.
        ymap[rho0==0] = 0.
        theta = self.theta
        ctheta = np.cos(theta)
        stheta = np.sin(theta)
        x = xmap*ctheta-ymap*stheta
        y = ymap*ctheta+xmap*stheta
        return x,y

        return self.align_coords(xmap,ymap,True)

    def alphax(self,x,y):
        pass

    def alphay(self,x,y):
        pass


class _PowerLaw(_MassModel):
    """
    A subclass for power-law mass models. The `power-law' aspect doesn't
        currently work, but it does work for SIE models.
    """
    def __init__(self,x=None,y=None,b=None,eta=1.,pa=None,q=None,load=False):
        self.x = x
        self.y = y
        self.b = b
        self.eta = eta
        self.pa = pa
        self.q = q

    def kappa(self, xin, yin):
        R = np.sqrt(self.q**2*xin*xin+yin*yin)
        rho = (self.b/R)**self.eta
        return 0.5*rho**(self.eta-2.)

    def deflections(self,xin,yin):
        if self.NoFreeParams==True:
            try:
                return self.deflx,self.defly
            except: pass

        x,y = self.align_coords(xin,yin)
        q = self.q
        if q==1.:
            q = 1.-1e-7  # Avoid divide-by-zero errors
        eps = (1.-q**2)**0.5
        if self.eta==1.:
            # SIE models
            r = (x**2+y**2)**0.5
            r[r==0.] = 1e-9
            xout = self.b*np.arcsinh(eps*x/q/r)*q**0.5/eps
            yout = self.b*np.arcsin(eps*y/r)*q**0.5/eps
        else:
            from powerlaw import powerlawdeflections as pld
            b,eta = self.b,self.eta
            s = 1e-7
            if x.ndim>1:
                yout,xout = pld(-1*y.ravel(),x.ravel(),b,eta,s,q)
                xout,yout = xout.reshape(x.shape),-1*yout.reshape(y.shape)
            else:
                yout,xout = pld(-1*y,x,b,eta,s,q)
                yout = -1*yout
        theta = -(self.theta - np.pi/2.)
        ctheta = np.cos(theta)
        stheta = np.sin(theta)
        x = xout*ctheta+yout*stheta
        y = yout*ctheta-xout*stheta
        self.deflx=x
        self.defly=y
        return x,y

    def shears(self,rho):
        return 0.0

    def phi(self, rho):
        return 0.0


class _ExtShear(_MassModel):

    def __init__(self,x=None,y=None,b=0.,theta=0.,pa=None):
        self.x = x
        self.y = y
        self.b = b
        if pa is None:
            self.theta = theta
        else:
            self.pa = pa

    def deflections(self,x,y):
        x = x-self.x
        y = y-self.y
        theta = self.theta #- pi/2
        s = np.sin(2*theta)
        c = np.cos(2*theta)

        # From Kormann B1422 paper
        alpha_x = -self.b*(x*c+y*s)
        alpha_y = -self.b*(x*s-y*c)

        return alpha_x,alpha_y


class _MassSheet(_MassModel):

    def __init__(self,x=None,y=None,b=0.):
        self.x = x
        self.y = y
        self.b = b

    def deflections(self,x,y):
        x = x-self.x
        y = y-self.y

        #\vec(alpha)=-kappa*\vec(x)
        alpha_x = -self.b*x
        alpha_y = -self.b*y

        return alpha_x,alpha_y


if __name__ == '__main__':
    mass_epl = _PowerLaw(x=0.0,y=0.0,b=1.5,eta=1.,pa=130.,q=0.7)
    print mass_epl.kappa(0.5,0.7)

