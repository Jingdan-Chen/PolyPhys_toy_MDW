import json
from typing import Dict, Any

def read_and_validate_config(file_path: str) -> Dict[str, Any]:
    """
    Read a polymerization configuration JSON file and validate the types.
    
    Args:
        file_path: Path to the JSON configuration file
        
    Returns:
        Dictionary containing the configuration parameters
        
    Raises:
        ValueError: If any parameter has an incorrect type
    """
    # Read the JSON file
    with open(file_path, 'r') as f:
        # Load JSON with comments (need to strip them first)
        file_content = ""
        for line in f:
            # Remove comments (everything after //)
            cleaned_line = line.split('//')[0]
            file_content += cleaned_line
        
        # Parse the cleaned JSON
        try:
            config = json.loads(file_content)
        except json.JSONDecodeError as e:
            raise ValueError(f"Invalid JSON format: {e}")
    
    # Define expected types for each parameter
    expected_types = {
        "I0": float,
        "n_max": int,
        "ending_thresh": float,
        "max_time": float,
        "dt": float,
        "flag1": bool,
        "flag2": bool,
        "flag3": bool,
        "flag4": bool,
        "flag5": bool,
        "x": int,
        "y": int,
        "T": float
    }
    
    # Validate parameter types
    for param, expected_type in expected_types.items():
        if param not in config:
            raise ValueError(f"Missing parameter: {param}")
        
        if not isinstance(config[param], expected_type):
            raise ValueError(f"Parameter '{param}' should be of type {expected_type.__name__}, but got {type(config[param]).__name__}")
        
        # require all float parameters to be positive
        if expected_type == float and config[param] < 0:
            raise ValueError(f"Parameter '{param}' should be positive/0")
    
    
    
    if config['n_max'] <= 0:
        raise ValueError("Parameter 'n_max' should be positive")
    
    if config["x"] not in [1, 2]:
        raise ValueError("Parameter 'x' should be either 1 or 2")
    
    if config["y"] not in [1, 2]:
        raise ValueError("Parameter 'y' should be either 1 or 2")
    
    return config