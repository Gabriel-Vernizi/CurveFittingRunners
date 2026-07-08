import pandas as pd
import scipy as scp
import numpy as np

def modelo_sigmoide(T, A, k, t0, V_base):
    return V_base + A / (1.0 + np.exp(k * (T - t0)))

def modelo_potencia(T, S, E):
    return S * T**(E - 1.0)

def modelo_hyper(T, a, alpha, gamma):
    # return (a * T**alpha) / (T**alpha + gamma**alpha) + a 
    return (a * gamma**alpha) / (T**alpha + gamma**alpha)