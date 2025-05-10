import math

T0 = 273.15  # K
T_exp = 50 + T0  # K
# R = 8.314  # J/(mol*K)
R = 1.987  # cal/(mol*K)
NA = 6.022E23  # mol^-1

def get_ki(T=T_exp):
    return 6.23E16 * math.exp(-3.0704E4/R/T)   # L/mol/min

def get_density(T=T_exp):
    return 0.9665 - 0.001164 * (T - T0) # g/cm^3

def get_kp(T=T_exp):
    return 2.952E7 * math.exp(-4343/R/T)  # L/mol/min

def get_kp_(dH=-143.6, T=T_exp): # cal/mol
    print(math.exp(dH/R/T))
    return get_kp(T) * math.exp(dH/R/T)  # L/mol/min

def get_kt(T=T_exp):
    return 5.88E9 * math.exp(-701/R/T)  # L/mol/min

def get_ktd_ratio(T=T_exp):
    return 1.6093 * math.exp(-440.12/R/T)  # L/mol/min

def get_ktd(T=T_exp):
    return get_kt(T) * get_ktd_ratio(T)  # L/mol/min

def get_ktc(T=T_exp):
    return get_kt(T) * (1 - get_ktd_ratio(T))  # L/mol/min

def get_kfm(T=T_exp):
    return 9.3435E4 * math.exp(-7475/R/T)  # L/mol/min

def get_f(T=T_exp):
    return 0.0247 * math.exp(2166 /R/T)  # L/mol/min

def get_kit(T=T_exp):
    # 2.26*10-6exp (-6578/RT)
    return 2.26E-6 * math.exp(-6578/R/T)  # L/mol/min

def get_Mw():
    return 100.12  # g/mol