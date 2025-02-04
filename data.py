import math

scenario = {
    "Fig. 1": {
        "title": "Matador, Saskatchewan",
        "year": 1971,
        "month": 7,
        "day": 30,
        "phi": 50.80 * (math.pi / 180),
        "albedo": 0.2,
        "lambda_e": -107.951 * (math.pi / 180),
        "sigma_C_H": lambda t :0,
        "sigma_C_M": lambda t: 0.06 if (17 < t < 22) else 0,
        "sigma_C_L": lambda t: 0.06 if (20 < t < 22) else 0,
        
        # Uncertainties :
        "D_j_t": 1,
        "D_j_r":1,
        "D_J":1,
        "Dt": 0.05,
        "D_sigma_C_H": 0.2,
        "D_sigma_C_M": 0.2,
        "D_sigma_C_L": 0.2,
        "D_S0": 20,
        "D_alpha":0.15
    },
        
        
        
    "Fig. 2": {
        "title": "Mizuho Station, Antarctica",
        "year": 1979,
        "month":11,
        "day": 13,
        "phi": -70.698* (math.pi / 180),
        "albedo": 0.8,
        "lambda_e": 44.332 * (math.pi / 180),
        "sigma_C_H": lambda t :0,
        "sigma_C_M": lambda t :0,
        "sigma_C_L": lambda t :0,
        
        # Uncertainties :
        "D_j_t": 1,
        "D_j_r":1,
        "D_J":1,
        "Dt": 0.05,
        "D_sigma_C_H": 0.2,
        "D_sigma_C_M": 0.2,
        "D_sigma_C_L": 0.2,
        "D_S0": 20,
        "D_alpha":0.15        
    },
        
        
        
    "Fig. 3": {
        "title": "Lake Ontario, Grimsby, Ontario",
        "year": 1969,
        "month": 8,
        "day": 28,
        "phi": 43.200 * (math.pi / 180),
        "albedo": 0.12,
        "lambda_e": -79.562 * (math.pi / 180),
        "sigma_C_H": lambda t :0,
        "sigma_C_M": lambda t :0,
        "sigma_C_L": lambda t :0,
                
        # Uncertainties :
        "D_j_t": 1,
        "D_j_r":1,
        "D_J":1,
        "Dt": 0.05,
        "D_sigma_C_H": 0.2,
        "D_sigma_C_M": 0.2,
        "D_sigma_C_L": 0.2,
        "D_S0": 20,
        "D_alpha":0.15
    },
        
        
        
    "Fig. 4": {
        "title": "Rothamsted, UK",
        "year": 1963,
        "month": 7,
        "day": 23,
        "phi": 51.809 * (math.pi / 180),
        "albedo": 0.22,
        "lambda_e": -0.355 * (math.pi / 180),
        "sigma_C_H": lambda t: 0.133 if (12 < t < 14) else 0,
        "sigma_C_M": lambda t: 0.5 if (12 < t < 14) else 0,
        "sigma_C_L": lambda t: 0.133 if (12 < t < 14) else 0,
                
        # Uncertainties :
        "D_j_t": 1,
        "D_j_r":1,
        "D_J":1,
        "Dt": 0.05,
        "D_sigma_C_H": 0.2,
        "D_sigma_C_M": 0.2,
        "D_sigma_C_L": 0.2,
        "D_S0": 20,
        "D_alpha":0.15
    },
        
        
        
    "Fig. 5": {
        "title": "Cedar River, Washington",
        "year": 1972,
        "month": 8,
        "day": 10,
        "phi": 47.47 * (math.pi / 180),
        "albedo": 0.15,
        "lambda_e": -122.16 * (math.pi / 180),
        "sigma_C_H": lambda t :0,
        "sigma_C_M": lambda t: 0.06 if (6 < t < 13) else 0,
        "sigma_C_L": lambda t: 0.12 if (6 < t < 10) else 0,
                
        # Uncertainties :
        "D_j_t": 1,
        "D_j_r":1,
        "D_J":1,
        "Dt": 0.05,
        "D_sigma_C_H": 0.2,
        "D_sigma_C_M": 0.2,
        "D_sigma_C_L": 0.2,
        "D_S0": 20,
        "D_alpha":0.15
    },
}

SIGMA = 5.67e-8  
EPS_CB = 0.85    
T_CB = 283.15   
T_SURFACE = 288.15  