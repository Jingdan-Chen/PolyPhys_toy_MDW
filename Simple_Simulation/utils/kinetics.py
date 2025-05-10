import numpy as np
import warnings
from multiprocessing import Pool

def decayed_coef(exponent=.05, rou=3., N=200):
    """
    Calculate a decay coefficient array for chain-length-dependent rate constants.
    
    This function creates a decay factor that accounts for the decrease in reactivity
    as the polymer chain grows longer. It follows an exponential decay model.
    
    Args:
        exponent: Power law exponent controlling decay rate (default: 0.05)
        rou: Coefficient controlling the magnitude of decay (default: 3.0)
        N: Length of the array to generate (default: 200)
        
    Returns:
        numpy.ndarray: Array of decay coefficients
    """
    array = np.linspace(1, N, N)
    res = np.exp(rou * (np.power(array, -exponent) - 1))
    return res
    # return np.ones(N)  # Uncomment to disable decay effect

def pad_array_to_length(arr: np.ndarray, n_shape: int, padding_last=True) -> np.ndarray:
    """
    Pad a 1D NumPy array to a specified length with zeros.
    
    Args:
        arr: A 1D NumPy array to be padded
        n_shape: The desired length of the padded array
        padding_last: If True, pad at the end; if False, pad at the beginning
        
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
    """
    Calculate the rate of change of initiator concentration.
    
    Args:
        It: Initiator concentration
        R_arr: Array of radical concentrations
        P_arr: Array of dead polymer concentrations
        M1: Monomer concentration
        config: Dictionary containing reaction parameters
        
    Returns:
        numpy.ndarray: Rate of change of initiator concentration
    """
    return np.array(-It * config['ki'])  # Rate of initiator consumption

def M_kinetic(It, R_arr, P_arr, M1, config):
    """
    Calculate the rate of change of monomer concentration.
    
    Models monomer consumption through various mechanisms:
    - Propagation with growing radicals
    - Initiator-induced initiation
    - Self-initiation (if enabled)
    - Depolymerization (if enabled, adds monomer)
    - Chain transfer to monomer (if enabled)
    
    Args:
        It: Initiator concentration
        R_arr: Array of radical concentrations
        P_arr: Array of dead polymer concentrations
        M1: Monomer concentration
        config: Dictionary containing reaction parameters
        
    Returns:
        float: Rate of change of monomer concentration
    """
    # Extract configuration flags
    flag1 = config['flag1']  # self-initiation
    flag2 = config['flag2']  # depolymerization
    flag3 = config['flag3']  # CT to monomer
    flag4 = config['flag4']  # kp decay
    
    # Extract rate constants
    kp = config['kp']        # Propagation rate constant
    ki = config['ki']        # Initiation rate constant
    f = config['f']          # Initiator efficiency
    ktd = config['ktd'] if not config['flag5'] else config['ktd'] + config['ktc']
    kit = config['kit'] if flag1 else 0  # Self-initiation rate constant
    kp_ = config['kp_'] if flag2 else 0  # Depolymerization rate constant
    kfm = config['kfm'] if flag3 else 0  # Chain transfer to monomer rate constant
    x = config['x'] if flag1 else 0      # Self-initiation order
    
    # Calculate monomer consumption through propagation
    if flag4:  # Chain-length-dependent propagation
        result = -kp * M1 * np.dot(R_arr, decayed_coef(N=len(R_arr))) - ki * M1 * It * f 
    else:      # Chain-length-independent propagation
        result = -kp * M1 * R_arr.sum() - ki * M1 * It * f
    
    # Add contributions from other mechanisms
    if flag1:  # Self-initiation
        result += - x * kit * np.power(M1, x) 
    if flag2:  # Depolymerization (adds monomer)
        result += kp_ * R_arr[1:].sum()
    if flag3:  # Chain transfer to monomer
        result += - kfm * R_arr.sum() * M1  
        
    return result

def R1_kinetic(It, R_arr, P_arr, M1, config):
    """
    Calculate the rate of change of primary radical (R1) concentration.
    
    Models primary radical dynamics through:
    - Formation from initiator
    - Consumption through propagation
    - Termination with other radicals
    - Self-initiation (if enabled)
    - Depolymerization (if enabled)
    - Chain transfer to monomer (if enabled)
    
    Args:
        It: Initiator concentration
        R_arr: Array of radical concentrations
        P_arr: Array of dead polymer concentrations
        M1: Monomer concentration
        config: Dictionary containing reaction parameters
        
    Returns:
        float: Rate of change of primary radical concentration
    """
    # Extract configuration flags
    flag1 = config['flag1']  # self-initiation
    flag2 = config['flag2']  # depolymerization
    flag3 = config['flag3']  # CT to monomer
    flag5 = config['flag5']  # ktc depreciation
    
    # Extract rate constants and parameters
    ki = config['ki']        # Initiation rate constant
    f = config['f']          # Initiator efficiency
    x = config['x'] if flag1 else 0  # Self-initiation order
    y = config['y'] if flag1 else 0  # Self-initiation stoichiometry
    kp = config['kp']        # Propagation rate constant
    kp_ = config['kp_'] if flag2 else 0  # Depolymerization rate constant
    ktc = config['ktc'] if not flag5 else 0  # Termination by combination
    ktd = config['ktd'] if not flag5 else config['ktd'] + config['ktc']  # Termination by disproportionation
    kfm = config['kfm'] if flag3 else 0  # Chain transfer to monomer
    kit = config['kit'] if flag1 else 0  # Self-initiation rate constant
    
    # Get current concentrations of key species
    R1 = R_arr[0]  # Primary radical
    R2 = R_arr[1]  # Secondary radical
    
    # Calculate R1 dynamics
    # Formation from initiator - consumption by propagation - termination with all radicals
    result = ki * M1 * It * f - kp * M1 * R1 - (ktc + ktd) * R1 * (R_arr.sum())
        
    # Add contributions from other mechanisms
    if flag1:  # Self-initiation
        result += y * np.power(M1, x) * kit
    if flag2:  # Depolymerization
        result += kp_ * R2    # R2 depolymerization produces R1
    if flag3:  # Chain transfer to monomer
        result += kfm * M1 * R_arr[1:].sum()  # CT creates primary radicals
    
    return result

def Rn_kinetic(It, R_arr, P_arr, M1, config):
    """
    Calculate the rate of change of secondary radical (Rn, n>1) concentrations.
    
    Models secondary radical dynamics through:
    - Formation through propagation of smaller radicals
    - Consumption through propagation to larger radicals
    - Termination with other radicals
    - Depolymerization (if enabled)
    - Chain transfer to monomer (if enabled)
    
    Args:
        It: Initiator concentration
        R_arr: Array of radical concentrations
        P_arr: Array of dead polymer concentrations
        M1: Monomer concentration
        config: Dictionary containing reaction parameters
        
    Returns:
        numpy.ndarray: Rates of change for each secondary radical concentration
    """
    # Extract configuration flags
    flag2 = config['flag2']  # depolymerization
    flag3 = config['flag3']  # CT to monomer
    flag4 = config['flag4']  # kp decay
    flag5 = config['flag5']  # ktc depreciation
    
    # Extract rate constants
    kp = config['kp']        # Propagation rate constant
    kp_ = config['kp_'] if flag2 else 0  # Depolymerization rate constant
    ktc = config['ktc'] if not flag5 else 0  # Termination by combination
    ktd = config['ktd'] if not flag5 else config['ktd'] + config['ktc']  # Termination by disproportionation
    kfm = config['kfm'] if flag3 else 0  # Chain transfer to monomer
    
    # Get radical concentrations for different chain lengths
    Rn = R_arr[1:]        # All secondary radicals (n>1)
    Rn_p1 = pad_array_to_length(R_arr[2:], Rn.shape[0])  # Radicals with length n+1
    Rn_m1 = R_arr[:-1]    # Radicals with length n-1
    
    # Calculate Rn dynamics with or without chain-length-dependent propagation
    if flag4:  # Chain-length-dependent propagation
        result = kp * M1 * Rn_m1 * decayed_coef(N=len(R_arr)-1) \
            - kp * M1 * Rn * decayed_coef(N=len(R_arr)+1)[2:] \
            - (ktc + ktd) * (Rn * R_arr.sum()) \
            - (ktc + ktd) * np.power(Rn, 2)
    else:      # Chain-length-independent propagation
        result = kp * M1 * Rn_m1 - kp * M1 * Rn \
            - (ktc + ktd) * (Rn * R_arr.sum()) \
            - (ktc + ktd) * np.power(Rn, 2)
        
    # Add contributions from other mechanisms
    if flag2:  # Depolymerization
        result += kp_ * (Rn_p1 - Rn)  # Rn+1 depropagates to Rn, Rn depropagates to Rn-1
    if flag3:  # Chain transfer to monomer
        result += - kfm * M1 * Rn     # Chain transfer consumes Rn
        
    return result

def zero_upper_diagonal(matrix):
    """
    Zero out the upper diagonal elements of a matrix.
    
    This helper function is used in termination by combination calculations.
    
    Args:
        matrix: Input matrix
        
    Returns:
        numpy.ndarray: Matrix with upper diagonal elements set to zero
    """
    result = matrix.copy()
    n, m = result.shape
    mask = np.triu(np.ones((n, m)), k=1)
    result[mask.astype(bool)] = 0
    
    return result

def calculate_trace(args):
    """
    Calculate the trace of a matrix with the given offset.
    
    This is a helper function for parallel processing of matrix traces.
    
    Args:
        args: Tuple containing (matrix, offset)
        
    Returns:
        float: Trace of the matrix with the given offset
    """
    mat, offset = args
    return np.trace(mat, offset=offset)

def Pn_kinetic(It, R_arr, P_arr, M1, config, thresh=1E-3):
    """
    Calculate the rate of change of dead polymer (Pn) concentrations.
    
    Models dead polymer formation through:
    - Termination by disproportionation
    - Termination by combination (if enabled)
    - Chain transfer to monomer (if enabled)
    
    Args:
        It: Initiator concentration
        R_arr: Array of radical concentrations
        P_arr: Array of dead polymer concentrations
        M1: Monomer concentration
        config: Dictionary containing reaction parameters
        thresh: Threshold for numerical stability (default: 1E-3)
        
    Returns:
        numpy.ndarray: Rates of change for each dead polymer concentration
    """
    # Extract configuration flags
    flag3 = config['flag3']  # CT to monomer
    flag5 = config['flag5']  # ktc depreciation
    
    # Extract rate constants
    ktd = config['ktd'] if not flag5 else config['ktd'] + config['ktc']  # Termination by disproportionation
    ktc = config['ktc']  # Termination by combination
    kfm = config['kfm'] if flag3 else 0  # Chain transfer to monomer
    
    # Calculate polymer formation by termination via disproportionation
    result_td = ktd * R_arr * (R_arr.sum()) + ktd * np.power(R_arr, 2)
    
    # Calculate polymer formation by termination via combination (if enabled)
    if not flag5:
        result_tc = np.zeros_like(result_td)
        
        # Create matrix of Ri*Rj products
        mat1 = np.matmul(R_arr.reshape(-1, 1), R_arr.reshape(1, -1))
        mat2 = zero_upper_diagonal(mat1)[:, ::-1]
        n = len(result_td)
        
        # Calculate the contribution to each chain length by combination
        for i in range(1, len(result_td)):
            diag_offset = len(result_td) - i
            result_tc[i] = np.trace(mat2, offset=diag_offset)
    else:
        result_tc = 0  # Combination terminated disabled
    
    # Combine contributions from different mechanisms
    result = result_tc * ktc + result_td
    
    # Add contribution from chain transfer to monomer (if enabled)
    if flag3:
        result += kfm * R_arr * M1
    
    return result
    
def evolution(It, M1, R_arr, P_arr, dt, t, config):
    """
    Evolve the polymerization system one time step using the kinetic equations.
    
    This is the main integration function that calculates the changes in all species
    concentrations and updates them for the next time step.
    
    Args:
        It: Initiator concentration
        M1: Monomer concentration
        R_arr: Array of radical concentrations
        P_arr: Array of dead polymer concentrations
        dt: Time step
        t: Current time
        config: Dictionary containing reaction parameters
        
    Returns:
        tuple: (new_It, new_M, new_R_arr, new_Pn, new_t) - Updated concentrations and time
    """
    # Calculate rates of change for all species
    d_It = I_kinetic(It, R_arr, P_arr, M1, config)
    d_M = M_kinetic(It, R_arr, P_arr, M1, config)
    d_R1 = R1_kinetic(It, R_arr, P_arr, M1, config)
    d_Rn = Rn_kinetic(It, R_arr, P_arr, M1, config)
    d_Pn = Pn_kinetic(It, R_arr, P_arr, M1, config)
    
    # Validate shapes for consistency
    assert d_It.shape == It.shape, f"Shape mismatch for d_It, {d_It.shape} != {It.shape}"
    assert d_Rn.shape == R_arr[1:].shape, f"Shape mismatch for d_Rn, {d_Rn.shape} != {R_arr[1:].shape}"
    assert d_Pn.shape == P_arr.shape, f"Shape mismatch for d_Pn, {d_Pn.shape} != {P_arr.shape}"
    assert len(d_R1.shape) == 1, f"Shape mismatch for d_R1, {d_R1.shape} != {R_arr[0].shape}"
    assert len(d_M.shape) == 1, f"Shape mismatch for d_M, {d_M.shape} != {M1.shape}"
    
    # Update all concentrations using explicit Euler method
    new_It = It + d_It * dt
    new_M = M1 + d_M * dt
    new_R1 = R_arr[0] + d_R1 * dt
    new_Rn = R_arr[1:] + d_Rn * dt
    new_Pn = P_arr + d_Pn * dt
    
    # Combine R1 and Rn into a single array
    new_R_arr = np.concatenate((new_R1.reshape(-1, 1), new_Rn.reshape(-1, 1)), axis=0).reshape(-1)
    
    # Ensure no negative concentrations (numerical stability)
    if (new_R_arr < 0).any():
        new_R_arr[new_R_arr < 0] = 0
        warnings.warn("Negative values in R_arr after evolution, setting to zero")
    
    # Return updated state
    return new_It, new_M, new_R_arr, new_Pn, t + dt
