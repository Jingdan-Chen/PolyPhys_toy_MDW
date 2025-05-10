import json
import numpy as np
import os

# Define the path to the JSON file
file_path = 'exam_json/big_bas4.json'
traj_dir_name = "result/bMMA4"

from utils import read_and_validate_config

# # Example usage
# if __name__ == "__main__":
try:
    config = read_and_validate_config(file_path)
    print("Configuration validated successfully!")
    print(config)
except ValueError as e:
    print(f"Error: {e}")
    
# recalc all params
from Sty_AIBN_param import get_ki, get_density, get_kp_, get_kp, get_kit
from Sty_AIBN_param import get_kp_, get_kt, get_kfm, get_f, get_Mw, get_ktd_ratio
from Sty_AIBN_param import get_Tgm, get_Tgp, get_Cpm, get_Cpp, get_density_p
from Sty_AIBN_param import get_delta, get_Vfc, get_alpha_m, get_alpha_p
from Sty_AIBN_param import get_A, get_n, get_m, get_K3, get_l0
from utils.kinetics import get_Vf

T0 = 273.15
T = T0 + config['T']
config['ki'] = get_ki(T=T)    
config['density'] = get_density(T=T)
config['kp'] = get_kp(T=T)
config['kp_'] = get_kp_(T=T)
config['kt0'] = get_kt(T=T)
config['ktd_ratio'] = get_ktd_ratio(T=T)

config['kfm'] = get_kfm(T=T)
config['f'] = get_f(T=T)
config['kit'] = get_kit(T=T)
config['Mw'] = get_Mw()
config['Tg_m'] = get_Tgm()
config['Tg_p'] = get_Tgp()
config['Cpm'] = get_Cpm()
config['Cpp'] = get_Cpp()
config['density_p'] = get_density_p(T=T)
config['delta'] = get_delta()
config['Vfc'] = get_Vfc(T=T)
config['alpha_m'] = get_alpha_m()
config['alpha_p'] = get_alpha_p()
config['A'] = get_A()
config['n'] = get_n()
config['m'] = get_m()
config['K3'] = get_K3(T=T)
config['l0'] = get_l0()

# for transport diffusion control
config['Mw_cr'] = None
config['kt_cr'] = None

# Correction for AIBN -> (AIBN1/2)2
config['I0'] = config['I0'] * 2.
config['ki'] = config['ki'] * .5
config['M0'] = config['density'] / config['Mw'] * 1000.

from utils import evolution
import numpy as np
from time import sleep
import pickle

It = np.array([config['I0']])
R_arr = np.zeros((config['n_max'],))
P_arr = np.zeros((config['n_max'],))
M0 = config['M0']
mass0 = M0

dt = 1E-10 # Warmup
dt_max = config['dt']
# dt_max = 4E-6
max_time = config['max_time']
t = 0
count_round = 0

report_time = 1000 # dump every x steps
warmup_time = 90

Traj_results = []


if not os.path.exists(traj_dir_name):
    os.makedirs(traj_dir_name)
    
from time import sleep
while t < max_time:
    It, M0, R_arr, P_arr, t, config = evolution(It, M0, R_arr, P_arr, dt, t, config=config)
    
    # print(count_round)
    if count_round % report_time == 0:
        avg_P_length = np.sum( np.arange(1, len(P_arr)) * P_arr[1:] / np.sum(P_arr[1:]))
        avg_R_length = np.sum( np.arange(1, len(R_arr)+1) * R_arr / np.sum(R_arr))
        mass_balance = (np.arange(1, len(P_arr)+1) * P_arr).sum() +\
            (np.arange(1, len(R_arr)+1) * R_arr).sum() + M0[0]
        Vf = get_Vf(M0, P_arr, config)[0]
        print(f"Time: {t:.2f}, dt:{dt}, It: {It}, M: {M0[0]:.2f},R:{R_arr.sum()}, avg_P_length: {avg_P_length:.2f}, mass balance: {(mass_balance / mass0):.2f}, Vf/Vfc: {Vf/config['Vfc']:.2f}, avg_R_length: {avg_R_length:.2f}")
        # Traj_results.append([t, It, M0, P_arr, R_arr])
        # Save the current state to a file
        curr_state = {
            't': t,
            'dt': dt,
            'It': It.tolist(),
            'M0': M0[0],
            'R_arr': R_arr.tolist(),
            'P_arr': P_arr.tolist(),
            'kt': config['kt'],
            'Vf_Vfc': Vf/config['Vfc']
        }
        
        state_file_path = os.path.join(traj_dir_name, f"state_{count_round}.pkl")
        with open(state_file_path, 'wb') as f:
            pickle.dump(curr_state, f)
        
    if count_round % warmup_time == 0:
        dt = min(dt * 1.01, dt_max)
    count_round += 1
    
    if M0 < config['M0'] * .1:
        print("Conversion > 90%")
        break
    # sleep(1)
    
