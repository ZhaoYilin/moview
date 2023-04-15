import numpy as np
from scipy.special import factorial2 as fact2

__all__ = ['PrimitiveGaussian','ContractedGaussian']

class PrimitiveGaussian(object):
    """Primitive Gaussian functions.
    
    Attributes
    ----------
    origin : 3-list
        Coordinate of the nuclei.

    n : int
        Principal quantum number.
    
    shell : 3-tuple
        Angular momentum.
    
    coef : float
        Coefficent of Primitive Gaussian function.

    exp : float
        Primitive Gaussian exponent.
    
    Properties
    ----------
    norm: float
        Normalization factor.
    
    Methods
    -------
    __init__(self,type,origin,n,shell,coef,exp)
        Initialize the instance.

    """
    def __init__(self,coefficient,origin,shell,exponent):
        """Initialize the instance.

        Parameters
        ----------
        origin : list
            Coordinate of the nuclei.

        n : int
            Principal quantum number.
    
        shell : list
            Angular momentum.

        coef : float
            Primitive Gaussian coefficients.
    
        exp : float
            Primitive Gaussian exponent.
        """
        self.coefficient = coefficient
        self.origin = origin
        self.shell = shell
        self.exponent  = exponent

    def __call__(self,x,y,z):
        X = x-self.origin[0]
        Y = y-self.origin[1]
        Z = z-self.origin[2]
        rr = X**2+Y**2+Z**2
        result = np.power(X,self.shell[0])*\
            np.power(Y,self.shell[1])*\
            np.power(Z,self.shell[2])*\
            np.exp(-self.exponent*rr)
        return result

    @property
    def norm(self):
        """Normalization factors. 
        
        Return
        ------
        norm : list
            Normalization factors
        """
        i,j,k = self.shell
        norm = np.sqrt(np.power(2,2*(i+j+k)+1.5)*
                        np.power(self.exponent,i+j+k+1.5)/
                        fact2(2*i-1)/fact2(2*j-1)/
                        fact2(2*k-1)/np.power(np.pi,1.5))
        return norm


class ContractedGaussian(object):
    """Predetermined linear combination of radial parts of GTOs
    Atomic orbtial represented by Contracted Gaussian functions.
    
    Attributes
    ----------
    origin : 3-list
        Coordinate of the nuclei.

    n : int
        Principal quantum number.
    
    shell : 3-tuple
        Angular momentum.
    
    coefs : list
        Primitive Gaussian coefficients.

    exps : list
        Primitive Gaussian exponents.
    
    Properties
    ----------
    norm: list
        Normalization factor.
    
    Methods
    -------
    __init__(self,type,origin,n,shell=(),coefs=[],exps=[])
        Initialize the instance.

    """
    def __init__(self,coefficients=[],origin=[],shell=[],exponents=[]):
        """Initialize the instance.

        Parameters
        ----------
        origin : list
            Coordinate of the nuclei.

        n : int
            Principal quantum number.
    
        shell : list
            Angular momentum.

        coefs : list
            Primitive Gaussian coefficients.
    
        exps : list
            Primitive Gaussian exponents.
        """
        self.coefficients = coefficients
        self.origin = origin
        self.shell = shell
        self.exponents  = exponents

    def __call__(self,x,y,z):
        X = x-self.origin[0]
        Y = y-self.origin[1]
        Z = z-self.origin[2]
        rr = X**2+Y**2+Z**2

        cg = np.zeros(x.shape)
        for coef, exp in zip(self.coefficients,self.exponents):
            cg += np.power(X,self.shell[0])*\
                np.power(Y,self.shell[1])*\
                np.power(Z,self.shell[2])*\
                np.exp(-exp*rr)
        return cg

    @property
    def norm(self):
        """Normalization factors. 
        
        Return
        ------
        norm : list
            Normalization factors
        """
        i,j,k = self.shell
        # self.norm is a list of length equal to number primitives
        norm = np.sqrt(np.power(2,2*(i+j+k)+1.5)*
                        np.power(self.exponents,i+j+k+1.5)/
                        fact2(2*i-1)/fact2(2*j-1)/
                        fact2(2*k-1)/np.power(np.pi,1.5))
        return norm

    @property
    def expansion(self):
        """Normalization factors. 
        
        Return
        ------
        primitives : list
            Normalization factors
        """
        n = len(self.coefficients)
        primitives = []
        for i in range(n):
            coefficient = self.coefficients[i]
            exponent = self.exponents[i]
            pg = PrimitiveGaussian(coefficient,self.origin,self.shell,exponent)
            primitives.append(pg)
        return primitives
