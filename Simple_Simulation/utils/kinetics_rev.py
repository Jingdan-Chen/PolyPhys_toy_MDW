import numpy as np
import warnings
from multiprocessing import Pool

def decayed_coef(exponent=.05, rou=3., N=200):
# def decayed_coef(exponent=2, rou=3., N=200):
    array = np.linspace(1, N, N)
    res = np.exp( rou * (np.power(array, -exponent) - 1) )
    return res
    # return np.ones(N)

def pad_array_to_length(arr: np.ndarray, n_shape: int, padding_last=True) -> np.ndarray:
    """
    Pad a 1D NumPy array to a specified length with zeros at the end.
    
    Args:
        arr: A 1D NumPy array to be padded
        n_shape: The desired length of the padded array
        
    Returns:
        A new NumPy array padded with zeros to the specified length
        
    Raises:
        ValueError: If the input array is not 1D or if n_shape is smaller than the length of the input array
    """
    # Validate input array is 1D
    if arr.ndim != 1:
        raise ValueError("Input array must be 1-dimensional")
    
    # Get the current length of the array
    current_length = arr.shape[0]
    
    # Check if padding is needed
    if current_length >= n_shape:
        return arr[:n_shape]  # Return array truncated or as is if n_shape <= current_length
    
    # Calculate the padding needed
    padding_needed = n_shape - current_length
    
    # Create a new array with zeros
    padded_array = np.zeros(n_shape, dtype=arr.dtype)
    
    # Copy the original array values
    if padding_last:
        padded_array[:current_length] = arr
    else:
        padded_array[1:] = arr
    
    return padded_array

def I_kinetic(It, R_arr, P_arr, M1, config):
    return np.array(-It * config['ki'])

def M_kinetic(It, R_arr, P_arr, M1, config):
    flag1 = config['flag1'] # self-initiation
    flag2 = config['flag2'] # depolymerization
    flag3 = config['flag3'] # CT to monomer
    flag4 = config['flag4'] # kp decay
    kp = config['kp']
    ki = config['ki']
    f = config['f']
    ktd = config['ktd'] if not config['flag5'] else config['ktd'] + config['ktc']
    kit = config['kit'] if flag1 else 0
    kp_ = config['kp_'] if flag2 else 0
    kfm = config['kfm'] if flag3 else 0
    x = config['x'] if flag1 else 0
    
    # M1 = P_arr[0]
    # R1 = R_arr[0]
    if flag4:
        result = -kp * M1 * np.dot(R_arr, decayed_coef(N=len(R_arr))) - ki * M1 * It * f 
    else:
        result = -kp * M1 * R_arr.sum() - ki * M1 * It * f
        # + ktd * R1 * (R_arr.sum()) \
        # + ktd * np.power(R1, 2) \
        
    if flag1:
        result += - x * kit * np.power(M1, x) 
    if flag2:
        result += kp_ * R_arr[1:].sum()
    if flag3:
        result += - kfm * R_arr.sum() * M1  
        
    return result

def R1_kinetic(It, R_arr, P_arr, M1, config):
    flag1 = config['flag1']
    flag2 = config['flag2']
    flag3 = config['flag3']
    flag5 = config['flag5']
    ki = config['ki']
    f = config['f']
    x = config['x'] if flag1 else 0
    y = config['y'] if flag1 else 0
    kp = config['kp']
    kp_ = config['kp_'] if flag2 else 0
    ktc = config['ktc'] if not flag5 else 0
    ktd = config['ktd'] if not flag5 else config['ktd'] + config['ktc']
    kfm = config['kfm'] if flag3 else 0
    kit = config['kit'] if flag1 else 0
    
    
    # M1 = P_arr[0]
    R1 = R_arr[0]
    R2 = R_arr[1]
    
    result = ki * M1 * It * f - kp * M1 * R1 \
        - (ktc + ktd) * R1 * (R_arr.sum())
        
    if flag1:
        result += y * np.power(M1,x) * kit
    if flag2:
        result += kp_ * R2
    if flag3:
        result += kfm * M1 * R_arr[1:].sum()
    return result

def Rn_kinetic(It, R_arr, P_arr, M1, config):
    flag2 = config['flag2']
    flag3 = config['flag3']
    flag4 = config['flag4']
    flag5 = config['flag5']
    kp = config['kp']
    kp_ = config['kp_'] if flag2 else 0
    ktc = config['ktc'] if not flag5 else 0
    ktd = config['ktd'] if not flag5 else config['ktd'] + config['ktc']
    kfm = config['kfm'] if flag3 else 0
    
    # M1 = P_arr[0]
    Rn = R_arr[1:]
    Rn_p1 = pad_array_to_length(R_arr[2:], Rn.shape[0])
    Rn_m1 = R_arr[:-1]
    
    if flag4:
        result = kp * M1 * Rn_m1 * decayed_coef(N=len(R_arr)-1) \
            - kp * M1 * Rn * decayed_coef(N=len(R_arr)+1)[2:] \
            - (ktc + ktd) * ( Rn * R_arr.sum() ) \
            - (ktc + ktd) * np.power(Rn, 2)
    else:
        result = kp * M1 * Rn_m1 - kp * M1 * Rn \
            - (ktc + ktd) * (Rn * R_arr.sum()) \
            - (ktc + ktd) * np.power(Rn, 2)
        # - (ktc + ktd) * (Rn * Rn.reshape(-1,1)).sum(axis=0) \
        
    if flag2:
        result += kp_ * (Rn_p1 - Rn)
    if flag3:
        result += - kfm * M1 * Rn
        
    return result

def zero_upper_diagonal(matrix):
    result = matrix.copy()
    n, m = result.shape
    mask = np.triu(np.ones((n, m)), k=1)
    result[mask.astype(bool)] = 0
    
    return result

def calculate_trace(args):
    mat, offset = args
    return np.trace(mat, offset=offset)

def Pn_kinetic(It, R_arr, P_arr, M1, config, thresh=1E-3):
    flag3 = config['flag3']
    flag5 = config['flag5']
    ktd = config['ktd'] if not flag5 else config['ktd'] + config['ktc']
    ktc = config['ktc']
    kfm = config['kfm'] if flag3 else 0
    
    # Rn = R_arr[1:]
    
    result_td = ktd * R_arr * (R_arr.sum()) \
        + ktd * np.power(R_arr, 2) \
    
    if not flag5:
        result_tc = np.zeros_like(result_td)
        
        mat1 = np.matmul(R_arr.reshape(-1, 1), R_arr.reshape(1, -1))
        mat2 = zero_upper_diagonal(mat1)[:, ::-1]
        n = len(result_td)
        for i in range(1, len(result_td)):
            diag_offset = len(result_td) - i
            result_tc[i] = np.trace(mat2, offset=diag_offset)
    else:
        result_tc = 0
        
        
    result = result_tc * ktc + result_td
    
    if flag3:
        result += kfm * R_arr * M1
    
    return result
    
    
def evolution(It, M1, R_arr, P_arr, dt, t, config):
    d_It = I_kinetic(It, R_arr, P_arr, M1,  config)
    d_M = M_kinetic(It, R_arr, P_arr, M1, config)
    d_R1 = R1_kinetic(It, R_arr, P_arr, M1, config)
    d_Rn = Rn_kinetic(It, R_arr, P_arr, M1, config)
    d_Pn = Pn_kinetic(It, R_arr, P_arr, M1, config)
    
    assert d_It.shape == It.shape, f"Shape mismatch for d_It, {d_It.shape} != {It.shape}"
    assert d_Rn.shape == R_arr[1:].shape, f"Shape mismatch for d_Rn, {d_Rn.shape} != {R_arr[1:].shape}"
    assert d_Pn.shape == P_arr.shape, f"Shape mismatch for d_Pn, {d_Pn.shape} != {P_arr.shape}"
    assert len(d_R1.shape) == 1, f"Shape mismatch for d_R1, {d_R1.shape} != {R_arr[0].shape}"
    assert len(d_M.shape) == 1, f"Shape mismatch for d_M, {d_M.shape} != {M1.shape}"
    
    # try:
    # Update the values
    new_It = It + d_It * dt
    new_M = M1 + d_M * dt
    new_R1 = R_arr[0] + d_R1 * dt
    new_Rn = R_arr[1:] + d_Rn * dt
    new_Pn = P_arr + d_Pn * dt
    # print(d_M + d_R1 + d_Rn.sum() + d_Pn.sum())
    # print(d_M, d_R1, d_Rn.sum(), d_Pn.sum())
    
    # new_P_arr = np.concatenate((new_M.reshape(-1, 1), new_Pn.reshape(-1, 1)), axis=0).reshape(-1)
    new_R_arr = np.concatenate((new_R1.reshape(-1, 1), new_Rn.reshape(-1, 1)), axis=0).reshape(-1)
        
    if (new_R_arr < 0).any():
        new_R_arr[new_R_arr < 0] = 0
        # raise ValueError("Negative values in R_arr after evolution, setting to zero")
        warnings.warn("Negative values in R_arr after evolution, setting to zero")
    # except Exception as e:
    #     # print(f"Error during evolution: {e}")
    #     dt=1E-10
    #     # Update the values
    #     new_It = It + d_It * dt
    #     new_M = P_arr[0] + d_M * dt
    #     new_R1 = R_arr[0] + d_R1 * dt
    #     new_Rn = R_arr[1:] + d_Rn * dt
    #     new_Pn = P_arr[1:] + d_Pn * dt
    #     # print(d_M + d_R1 + d_Rn.sum() + d_Pn.sum())
    #     # print(d_M, d_R1, d_Rn.sum(), d_Pn.sum())
        
    #     new_P_arr = np.concatenate((new_M.reshape(-1, 1), new_Pn.reshape(-1, 1)), axis=0).reshape(-1)
    #     new_R_arr = np.concatenate((new_R1.reshape(-1, 1), new_Rn.reshape(-1, 1)), axis=0).reshape(-1)
        
    #     if (new_R_arr < 0).any():
    #         new_R_arr[new_R_arr < 0] = 0
    #         raise ValueError("Negative values in R_arr after evolution, setting to zero")
    #         # warnings.warn("Negative values in R_arr after evolution, setting to zero")
        
        
    return new_It, new_M, new_R_arr, new_Pn, t + dt
    
    
     