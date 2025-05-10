import math

T0 = 273.15  # K
T_exp = 50 + T0  # K
# R = 8.314  # J/(mol*K)
R = 1.987  # cal/(mol*K)
NA = 6.022E23  # mol^-1

def get_ki(T=T_exp):
    # 6.33*10^{16}exp(-3.0719*10^4/RT)
    return 6.33E16 * math.exp(-3.0719E4/R/T)  # L/mol/min

def get_density(T=T_exp):
    # 0.924-0.000918 (T-273.15)
    return 0.924 - 0.000918 * (T - T0)  # g/cm^3

def get_kp(T=T_exp):
    # 1.302*10^9 exp((-7759.2/RT))
    return 1.302E9 * math.exp(-7759.2/R/T)  # L/mol/min

def get_kp_(dH=-17000, T=T_exp): # cal/mol
    print(math.exp(dH/R/T))
    return get_kp(T) * math.exp(dH/R/T)  # L/mol/min

def get_kt(T=T_exp):
    # 4.92*10^{11} exp((-3471.3/RT)^1)
    return 4.92E11 * math.exp(-3471.3/R/T)  # L/mol/min

def get_ktd_ratio(T=T_exp):
    return 0.01

def get_ktd(T=T_exp):
    return get_kt(T) * get_ktd_ratio(T)  # L/mol/min

def get_ktc(T=T_exp):
    return get_kt(T) * (1 - get_ktd_ratio(T))  # L/mol/min

def get_kfm(T=T_exp):
    # 1.386*10^8 exp((-12670/RT))
    return 1.386E8 * math.exp(-12670/R/T)  # L/mol/min

def get_f(T=T_exp):
    return 0.6

def get_kit(T=T_exp):
    # 1.35*10^7 exp((-27450/RT))
    return 1.35E7 * math.exp(-27450/R/T)  # L/mol/min

def get_Mw():
    return 104.12  # g/mol