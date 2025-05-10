import math

T0 = 273.15  # K
T_exp = 50 + T0  # K
# R = 8.314  # J/(mol*K)
R = 1.987  # cal/(mol*K)
NA = 6.022E23  # mol^-1

def get_ki(T=T_exp):
    return 6.23E16 * math.exp(-3.0704E4/R/T)   # L/mol/min

def get_density(T=T_exp):
    # 0.919-0.001012\(T-273.15)
    return 0.919 - 0.001012 * (T_exp - T0)  # g/cm^3

def get_kp(T=T_exp):
    # 1.326*10^9exp(-4278.1/RT)
    return 1.326E9 * math.exp(-4278.1/R/T)  # L/mol/min

def get_kp_(dH=-18400, T=T_exp): # cal/mol
    print(math.exp(dH/R/T))
    return get_kp(T) * math.exp(dH/R/T)  # L/mol/min

def get_kt(T=T_exp):
    # 8.04*10^{10}exp(-1338.4/RT)^1
    return 8.04E10 * math.exp(-1338.4/R/T)  # L/mol/min

def get_ktd_ratio(T=T_exp):
    return .1

def get_ktd(T=T_exp):
    return get_kt(T) * get_ktd_ratio(T)  # L/mol/min

def get_ktc(T=T_exp):
    return get_kt(T) * (1 - get_ktd_ratio(T))  # L/mol/min

def get_kfm(T=T_exp):
    # 9.3436*10^5exp(-7475/RT)
    return 9.3436E5 * math.exp(-7475/R/T)  # L/mol/min

def get_f(T=T_exp):
    return 0.0247 * math.exp(2166 /R/T)  # L/mol/min

def get_kit(T=T_exp):
    # 4.96*10^4exp(-17483/RT)
    return 4.96E4 * math.exp(-17483/R/T)  # L/mol/min

def get_Mw():
    # 1.0E6
    return 128.17  # g/mol