from abc import ABC, abstractmethod
import ROOT as R
import numpy as np

class PDF(ABC):

    @abstractmethod
    def __init__(self) -> None:
        """
        
        """
        pass


    @abstractmethod
    def evaluate(self, x: float|np.ndarray, norm: float) -> float:
        """
        
        """
        pass


    @abstractmethod
    def get_integral(self, min: float, max: float) -> float:
        """
        
        """
        pass


    def get_normalisation(self, nbins: int, N: float, min: float, max: float) -> float:
        integral = self.get_integral(min, max)
        bin_width = (max - min) / nbins
        normalisation = bin_width * N / integral

        return normalisation



class GaussianPDF(PDF):

    def __init__(self, mu: float, sigma: float) -> None:
        self._mu = mu
        self._sigma = sigma


    def evaluate(self, x: float|np.ndarray, norm: float=1.) -> float|np.ndarray:
        A = norm * 1 / (self._sigma * np.sqrt(2*np.pi))
        G = np.exp(-0.5 * np.power((x - self._mu) / self._sigma, 2))

        return A * G


    def get_integral(self, min: float, max: float) -> float:
        x = R.RooRealVar("x", "x", min, max)
        mu = R.RooRealVar("mu", "mu", self._mu)
        sigma = R.RooRealVar("sigma", "sigma", self._sigma)
        G = R.RooGaussian("G", "G", x, mu, sigma)

        # must multiply by this because of the way
        # RooGaussians are coded in RooFit
        A = 1 / (self._sigma * np.sqrt(2*np.pi))

        return A * G.analyticalIntegral(1, "")
        


class BifurGaussPDF(PDF):

    def __init__(self, mu: float, sigmaL: float, sigmaR: float) -> None:
        self._mu = mu
        self._sigmaL = sigmaL
        self._sigmaR = sigmaR

        self._GaussianL = GaussianPDF(mu, sigmaL)
        self._GaussianR = GaussianPDF(mu, sigmaR)


    def evaluate(self, x: float|np.ndarray, norm: float=1.) -> float|np.ndarray:
        # must multiply out the normalisation of regular
        # GaussianPDFs because L and R have different widths
        L = (x < self._mu)  * (self._sigmaL * np.sqrt(2*np.pi)) * self._GaussianL.evaluate(x, norm)
        R = (x >= self._mu) * (self._sigmaR * np.sqrt(2*np.pi)) * self._GaussianR.evaluate(x, norm)
        
        return L + R


    def get_integral(self, min: float, max: float) -> float:
        x = R.RooRealVar("x", "x", min, max)
        mu = R.RooRealVar("mu", "mu", self._mu)
        sigmaL = R.RooRealVar("sigmaL", "sigmaL", self._sigmaL)
        sigmaR = R.RooRealVar("sigmaR", "sigmaR", self._sigmaR)
        BG = R.RooBifurGauss("BG", "BG", x, mu, sigmaL, sigmaR)

        return BG.analyticalIntegral(1, "")



class JohnsonPDF(PDF):

    def __init__(self, mu: float, sigma: float, gamma: float, delta: float) -> None:
        self._mu = mu
        self._sigma = sigma
        self._gamma = gamma
        self._delta = delta


    def evaluate(self, x: float|np.ndarray, norm: float=1.) -> float|np.ndarray:
        A = norm * self._delta / (self._sigma * np.sqrt(2*np.pi))
        z = self._gamma + self._delta * np.arcsinh((x - self._mu) / self._sigma)
        J = 1 / (np.sqrt(1 + ((x - self._mu) / self._sigma)**2)) * np.exp(-0.5*z**2)

        return A * J


    def get_integral(self, min: float, max: float) -> float:
        x = R.RooRealVar("x", "x", min, max)
        mu = R.RooRealVar("mu", "mu", self._mu)
        sigma = R.RooRealVar("sigma", "sigma", self._sigma)
        gamma = R.RooRealVar("gamma", "gamma", self._gamma)
        delta = R.RooRealVar("delta", "delta", self._delta)
        J = R.RooJohnson("J", "J", x, mu, sigma, gamma, delta)

        return J.analyticalIntegral(1, "")



class CrystalBallPDF(PDF):

    def __init__(self, mu: float, sigmaL: float, sigmaR: float, aL: float, nL: float, aR: float, nR: float) -> None:
        self._mu = mu
        self._sigmaL = abs(sigmaL)
        self._sigmaR = abs(sigmaR)
        self._aL = abs(aL)
        self._nL = nL
        self._aR = abs(aR)
        self._nR = nR


    # def evaluate(self, x: float|np.ndarray, norm: float=1.) -> float|np.ndarray:
    #     AL = (self._nL / abs(self._aL))**self._nL * np.exp(-0.5 * self._aL**2)
    #     AR = (self._nR / abs(self._aR))**self._nR * np.exp(-0.5 * self._aR**2)

    #     BL = self._nL / abs(self._aL) - abs(self._aL)
    #     BR = self._nR / abs(self._aR) - abs(self._aR)

    #     One   =  ((x - self._mu) / self._sigmaL <  -self._aL) * \
    #         AL * (BL - (x - self._mu) / self._sigmaL)**-self._nL
    #     Two   = (((x - self._mu) / self._sigmaL >= -self._aL) & ((x - self._mu) / self._sigmaL <= 0)) * \
    #         np.exp(-0.5 * ((x - self._mu) / self._sigmaL)**2)
    #     Three = (((x - self._mu) / self._sigmaL > 0) & ((x - self._mu) / self._sigmaR <= self._aR)) * \
    #         np.exp(-0.5 * ((x - self._mu) / self._sigmaR)**2)
    #     Four  =  ((x - self._mu) / self._sigmaR > self._aR) * \
    #         AR * (BR + (x - self._mu) / self._sigmaR)**-self._nR

    #     return norm * (One + Two + Three + Four)

    def evaluate(self, x: float|np.ndarray, norm: float=1.) -> float|np.ndarray:
        # if type(x) == float:
        if len(np.shape(x)) == 0:
            return norm * self.evaluate_point(float(x))
        else:
            x = np.array(x)
            new_x = np.zeros_like(x)
            for i in range(x.shape[0]):
                new_x[i] = norm * self.evaluate_point(x[i])
            return new_x
    def evaluate_point(self, x: float) -> float:
        t = x - self._mu
        t /= self._sigmaL if x < self._mu else self._sigmaR

        if t < -self._aL:
            return self.evaluate_tail(t, self._aL, self._nL)
        elif t <= self._aR:
            return np.exp(-0.5 * t**2)
        else:
            return self.evaluate_tail(-t, self._aR, self._nR)
    def evaluate_tail(self, t: float, a: float, n: float) -> float:
        A = (n / a)**n * np.exp(-0.5 * a**2)
        B = n / a - a
        return A * np.power(B - t,-n)



    def get_integral(self, min: float, max: float) -> float:
        x = R.RooRealVar("x", "x", min, max)
        mu = R.RooRealVar("mu", "mu", self._mu)
        sigmaL = R.RooRealVar("sigmaL", "sigmaL", self._sigmaL)
        sigmaR = R.RooRealVar("sigmaR", "sigmaR", self._sigmaR)
        aL = R.RooRealVar("aL", "aL", self._aL)
        nL = R.RooRealVar("nL", "nL", self._nL)
        aR = R.RooRealVar("aR", "aR", self._aR)
        nR = R.RooRealVar("nR", "nR", self._nR)
        CB = R.RooCrystalBall("CB", "CB", x, mu, sigmaL, sigmaR, aL, nL, aR, nR)

        return CB.analyticalIntegral(1, "")



class ArgusPDF(PDF):

    def __init__(self, c: float, p: float, x0: float=139.57061) -> None:
        self._c = c
        self._p = p
        self._x0 = x0


    def evaluate(self, x: float|np.ndarray, norm: float=1.) -> float|np.ndarray:
        # check x value(s) are physical
        t = (2 * self._x0 - x) / self._x0
        if type(x) == np.ndarray:
            assert (t < 1).all()
        else:
            assert t < 1

        A = norm * (2 * self._x0 - x)
        u = 1 - t**2
        B = u**self._p
        C = np.exp(self._c * u)

        return A * B * C
    

    def get_integral(self, min: float, max: float, n_steps: int=1_000) -> float:
        # RooFits analyticalIntegral for RooArgusBG
        # misbehaves so do it numerically 
        X = np.linspace(min, max, n_steps+1)
        dX = X[1] - X[0]

        integral = 0.
        for i in range(n_steps):
            integral += dX * (self.evaluate(X[i]) + self.evaluate(X[i+1])) / 2

        return integral



class ExponentialPDF(PDF):

    def __init__(self, a: float) -> None:
        self._a = a


    def evaluate(self, x: float|np.ndarray, norm: float=1.) -> float|np.ndarray:
        E = norm * np.exp(self._a * x)

        return E


    def get_integral(self, min: float, max: float) -> float:
        x = R.RooRealVar("x", "x", min, max)
        a = R.RooRealVar("a", "a", self._a)
        E = R.RooExponential("E", "E", x, a)

        return E.analyticalIntegral(1, "")
