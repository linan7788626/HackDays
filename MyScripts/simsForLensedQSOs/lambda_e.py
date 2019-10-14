import numpy as np
import pylab as pl
import scipy.special as ss
import scipy.integrate as sint


F = ss.hyp2f1
#--------------------------------------------------------------------
def Intergral_ud(llim,ulim,f,nbins=100000):

    soften = 1e-6
    #--------------------------------------
    # Scipy
    #
    # res = sint.romberg(f, llim+soften, ulim-soften)
    res = sint.quad(f, llim+soften, ulim-soften)[0]

    # #--------------------------------------
    # # Simplest
    # #
    # x = np.linspace(llim+soften, ulim-soften, nbins)
    # xl = x[:-1]
    # xu = x[1:]
    # dx = x[1] - x[0]

    # res = np.sum((f(xl)+f(xu))*dx/2.0)

    # #--------------------------------------
    # # Simpson's Rule
    # #
    # x = np.linspace(llim+soften, ulim-soften, nbins)
    # xl = x[:-1]
    # xu = x[1:]
    # xm = 0.5*(x[:-1]+x[1:])
    # dh = (x[1] - x[0])/2.0

    # res = np.sum((f(xl)+4.0*f(xm)+f(xu))*dh/3.0)

    return res


def lambda_e_obl(e):
    f2 = 1.0-e

    def WX_obl_qt(t, f2):
        # f2_p = np.sqrt(1.0-f2*f2)
        # q3 = np.sqrt(1.0-f2_p**2/t**2)
        # q3_p = np.sqrt(1.0-q3*q3)
        q3 = np.sqrt(1.0-(1.0-f2*f2)/t**2)
        q3_p = np.sqrt((1.0-f2**2)/t**2)
        # q3_p = np.sqrt(np.abs(1.0-q3*q3))

        res = q3/q3_p*np.arcsin(q3_p)
        return res

    def WZ_obl_qt(t, f2):

        # f2_p = np.sqrt(1.0-f2*f2)
        # q3 = np.sqrt(1.0-f2_p**2/t**2)
        # q3_p = np.sqrt(1.0-q3*q3)
        q3 = np.sqrt(1.0-(1.0-f2*f2)/t**2)
        q3_p = np.sqrt((1.0-f2**2)/t**2)
        # q3_p = np.sqrt(np.abs(1.0-q3*q3))

        # res = 0.5*q3/(q3_p*np.arctan2(q3_p,q3)) \
            # * Intergral_ud(0.0, np.pi,
                       # lambda theta: (1.0/np.sin(theta)
                                      # *(np.arctan2(q3_p,q3)**2.0-np.arctan2(q3_p*np.cos(theta),q3)**2.0)))

        res = 0.5*q3/(q3_p*np.arctan(q3_p/q3)) \
            * Intergral_ud(0.0, np.pi,
                       lambda theta: (1.0/np.sin(theta)
                                      *(np.arctan(q3_p/q3)**2.0
                                        -np.arctan(q3_p*np.cos(theta)/q3)**2.0)))

        return res

    res = np.pi/4.0*f2**1.5*F(1.0, 0.5, 2.0, f2**2) \
        / Intergral_ud(np.sqrt(1.0-f2*f2), 1.0,
                       lambda t: ((WX_obl_qt(t, f2)-WZ_obl_qt(t, f2))*t*t + WZ_obl_qt(t, f2))/np.sqrt(1.0-t*t))

    return res


def lambda_e_prol(e):
    f2 = 1.0-e

    def WX_prol_qt(t, f2):
        # f2_p = np.sqrt(1.0-f2*f2)
        # q3 = f2_p*f2_p/f2/f2/t/t+1.0
        # q3_p = np.sqrt(q3*q3-1.0)
        # q_p = |1 - q^2|^1/2
        q3 = np.sqrt(1.0+(1.0/f2**2-1.0)/t**2)
        q3_p = np.sqrt((1.0/f2**2-1.0)/t**2)

        res = 0.5*q3/q3_p*np.log((q3+q3_p)/(q3-q3_p))

        return res

    def WZ_prol_qt(t, f2):
        # f2_p = np.sqrt(1.0-f2*f2)
        # q3 = f2_p*f2_p/f2/f2/t/t+1.0
        # q3_p = np.sqrt(q3*q3-1.0)
        q3 = np.sqrt(1.0+(1.0/f2**2-1.0)/t**2)
        q3_p = np.sqrt((1.0/f2**2-1.0)/t**2)

        res = 0.25*q3/(q3_p*np.log((q3+q3_p)/(q3-q3_p))) \
        * Intergral_ud(0.0,np.pi,
                       lambda theta: (1.0/np.sin(theta) \
                                      * (np.log((q3+q3_p)/(q3-q3_p))**2 \
                                         - np.log((q3+q3_p*np.cos(theta)) \
                                                  /(q3-q3_p*np.cos(theta)))**2)))

        return res

    f2_p = np.sqrt(1.0-f2*f2)
    q3_max = 3.46717
    t_min = f2_p/f2/np.sqrt(q3_max*q3_max-1.0)

    res = Intergral_ud(t_min,1.0, \
                       lambda t: (np.sqrt((t*t*f2*f2+f2_p*f2_p) \
                                          /(f2*(1-t*t)))/t)) \
        / Intergral_ud(t_min,1.0, \
                       lambda t: (((WX_prol_qt(t, f2)-WZ_prol_qt(t, f2))*t*t \
                                   +WZ_prol_qt(t, f2))/np.sqrt(1-t*t)))

    return res


def lambda_e_tot(e, ef=0.5):
    res = ef*lambda_e_prol(e) + (1.0 - ef)*lambda_e_obl(e)
    return res


if __name__ == '__main__':
    et = np.linspace(1e-6, 0.6, 20)

    print lambda_e_obl(1e-10)
    print lambda_e_obl(1e-10)**2.0-1.0

    lto = np.array(map(lambda_e_obl, et))**2.0 - 1.0
    ltp = np.array(map(lambda_e_prol, et))**2.0 - 1.0
    lta = np.array(map(lambda_e_tot, et))**2.0 - 1.0

    pl.figure()
    pl.plot(et, lto, '-')
    pl.plot(et, ltp, '-')
    pl.plot(et, lta, '-')
    pl.show()
