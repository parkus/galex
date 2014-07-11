"""The variability module is a set of analysis tools for quantifying and
filtering the variability in time-series (or other sequential) data.

Modification history:
2014/05/12 - created

@author: Parke Loyd
"""
from numpy import sqrt, array, var, mean, exp, log, pi, inf
from scipy.optimize import minimize
from scipy.integrate import quad
import pdfutils as pu
import pdb
    
def excess_noise_PDF(y,base_noise,Poisson=False):
    """Computes the maximum likelihood value of the excess noise.
    
    Specifically, this function assumes the data, y, are each drawn from a
    normal distribution with the same mean but different variances. The
    variance for the PDF from which a given point was drawn is given by
    base_noise that may differ from point to point and some excess noise that
    is the same for all points.
    
    Modification History:
    2014/05/12 - Created (only for Poisson=False case)
    """
    
    y, base_noise = array(y), array(base_noise)
    ymn = mean(y)
    yvar = var(y)
    base_var = base_noise**2
    x_guess = sqrt(yvar - mean(base_var)) if yvar > mean(base_var) else 0.0
    
    def log_like(mn,x_noise,log_norm_fac):
        x_var = x_noise**2
        var = x_var + base_var
        terms = -log(2*pi*var)/2 - (y - mn)**2/2/var
        return sum(terms) + log_norm_fac
    
    #initial pass at a normalization factor - use peak value
    def neg_like(x):
        mn,x_noise= x
        return -log_like(mn,x_noise,0.0)
    result = minimize(neg_like, [ymn,x_guess], method='Nelder-Mead')
    log_norm_fac = result.fun
    mn_pk, x_noise_pk = result.x
    
    def like(mn,x_noise):
        return exp(log_like(mn, x_noise, log_norm_fac))
    
    def xmn_like(xnorm):
        constrained_like = lambda mn: like(mn, xnorm*mn)
        result = quad(constrained_like, mn_pk, inf)[0]
        result += -quad(constrained_like, mn_pk, 0.0)[0]
        return result
    
    xmn_pk = x_noise_pk/mn_pk
    result = minimize(lambda x: -xmn_like(x), xmn_pk, method='Nelder-Mead')
    xmn_pk = result.x[0]
#    f_pk = -result.fun
    if xmn_pk < 0.0: xmn_pk = 0.0
#    width = 1.0/f_pk
#    area = pu.gauss_integrate(xmn_like, xmn_pk, 0.0, inf, sigma_guess=width)
    area = quad(xmn_like, xmn_pk, inf)[0]
    area += -quad(xmn_like, xmn_pk, 0.0)[0]
    pdf = lambda x: xmn_like(x)/area if x >= 0.0 else 0.0
#    if pdf(xmn_pk) <= 0.0: pdb.set_trace()
    return pdf, xmn_pk