from random import *
import math
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.ndimage import gaussian_filter
from scipy.signal import savgol_filter
from data import *

def UTC(lamb):
    return lamb/math.pi*12

def t_shift(t, utc):
    """
    Shift time t by UTC.
    """
    if (t-utc) < 24 and (t-utc)>=0 :
        ts = t - utc
    elif (t - utc) >= 24:
        ts = (t - utc) - 24
    else:
        ts = 24 + (t-utc)
    return ts

def convert_decimal_hours_to_hms(decimal_hours):
    """
    Convert decimal hours to hours, minutes, and seconds.
    
    :param decimal_hours: Time in decimal hours.
    :return: Tuple containing (hours, minutes, seconds).
    """
    # Calculate hours
    hours = int(decimal_hours)
    # Calculate the remaining minutes
    remaining_minutes = (decimal_hours - hours) * 60
    minutes = int(remaining_minutes)
    # Calculate the remaining seconds
    seconds = (remaining_minutes - minutes) * 60
    return hours, minutes, int(seconds)

def time_of_year_from_date(date):
    """
    Calculate the time of the year j_t based on a UTC date.
    
    :param date: Datetime object representing the UTC date and time.
    :return: The day of the year as a decimal.
    """
    # Get the day of the year (from 1 to 365/366)
    day_of_year = date.timetuple().tm_yday
    
    # Calculate the decimal hour
    hour = date.hour 
    minute = date.minute
    second = date.second
    
    # Calculate time of the year including hours, minutes, and seconds
    j_t = day_of_year + (hour / 24) + (minute / (24 * 60)) + (second / (24 * 3600))
    return j_t


def solar_declination_angle(year, month, day, hour=0, minute=0, second=0):
    """
    Calculate the solar declination angle.
    """
    phi_r = 0.409 
    j_r = 173  
    j_t = time_of_year_from_date(datetime(year, month, day, hour, minute, second)) 
    J = 365  
    
    delta_s = phi_r * math.cos(2 * math.pi * (j_t - j_r) / J)
    return delta_s

def solar_elevation_angle(phi, delta_s, t, lambda_e):
    """
    param phi: latitude of the location
    param delta_s: solar declination angle
    param t: time of the day in hours
    param lambda_e: longitude of the location
    """
    sin_psi = math.sin(phi) * math.sin(delta_s) - math.cos(phi) * math.cos(delta_s) * math.cos((math.pi * t / 12) + lambda_e)
    return math.asin(sin_psi)


def incoming_shortwave_radiation(psi, sigma_CH, sigma_CM, sigma_CL):
    """
    :param psi: solar elevation angle (angle of sun above horizon)
    :param sigma_CH: cloud cover (0-1) for high level clouds
    :param sigma_CM: cloud cover (0-1) for middle level clouds
    :param sigma_CL: cloud cover (0-1) for low level clouds
    """
    S0 = 1370 # W/m^2 (solar constant)
    tau_K = (0.6 + 0.2 * math.sin(psi)) * (1 - 0.4 * sigma_CH) * (1 - 0.7 * sigma_CM) * (1 - 0.4 * sigma_CL)
    K_down = max(S0 * tau_K * math.sin(psi), 0)
    return K_down

def outgoing_shortwave_radiation(K_down, albedo):
    K_up = - albedo * K_down
    return K_up

def net_longwave_radiation(sigma_CH, sigma_CM, sigma_CL, mode = "BG"):
    if mode == "BG":
        L_star = -97.28 * (1 - 0.1 * sigma_CH - 0.3 * sigma_CM - 0.6 * sigma_CL)
    elif mode == "SPP":
        sigma_c = 1 - (1 - sigma_CH) * (1 - sigma_CM) * (1 - sigma_CL)
        correction = 1 - 0.94e-5 * T_SURFACE**2
        L_star = (
            -SIGMA * T_SURFACE**4 * correction
            + 0.3 * EPS_CB * SIGMA * T_CB**4 * sigma_c
        )
    return L_star

def radiation_budget(K_down, K_up, L_star):
    """net radiative flux"""
    R_N = K_down + K_up + L_star
    return R_N

def hourly_wawe_radiation(coordinates, albedo, date, cloud_cover):
    return incoming_shortwave_radiation, outgoing_shortwave_radiation, net_longwave_radiation, net_radiative_flux


def radiation_model(data,time):
    K_down_values,  K_up_values, L_star_values,  R_N_values = [], [], [], []
    utc = UTC(data["lambda_e"])
    time_shift = [t_shift(t, utc) for t in time]
    sigma_C_H = savgol_filter(gaussian_filter([data["sigma_C_H"](t) for t in time], 7), 51, 5) 
    sigma_C_M = savgol_filter(gaussian_filter([data["sigma_C_M"](t) for t in time], 7), 51, 5) 
    sigma_C_L = savgol_filter(gaussian_filter([data["sigma_C_L"](t) for t in time], 7), 51, 5)

    # plt.figure(figsize=(10, 6))
    # plt.plot(time, sigma_C_H, label='sigma_C_H')
    # plt.plot(time, sigma_C_M, label='sigma_C_M')
    # plt.plot(time, sigma_C_L, label='sigma_C_L')
    
    for i, t in enumerate(time_shift):

        # Convert time from decimal hours to hours, minutes, seconds
        h, m, s = convert_decimal_hours_to_hms(t)

        delta_s = solar_declination_angle(data["year"], data["month"], data["day"], h, m, s)

        psi = solar_elevation_angle(data["phi"], delta_s, t, data["lambda_e"])

        K_down = incoming_shortwave_radiation(psi, sigma_C_H[i], sigma_C_M[i], sigma_C_L[i])
        K_up = outgoing_shortwave_radiation(K_down, data["albedo"])
        L_star = net_longwave_radiation(sigma_C_H[i], sigma_C_M[i], sigma_C_L[i], mode = "BG")
        R_N = radiation_budget(K_down, K_up, L_star)

        K_down_values.append(K_down)
        K_up_values.append(K_up)
        L_star_values.append(L_star)
        R_N_values.append(R_N)

    return K_down_values, K_up_values, L_star_values, R_N_values, time_shift



#======================= Main program part ================================

time = np.arange(0, 24, 0.05)  # Time in hours

file_path = "radiation_profiles.xlsx"  
excel_data = pd.ExcelFile(file_path)

for sheet_name in excel_data.sheet_names:
    df = pd.read_excel(file_path, sheet_name=sheet_name)
    
    t = df["t "].iloc[1:-1].to_numpy()    
    k_down = df["K-down"].iloc[1:-1].to_numpy()  
    k_up = df["K-up"].iloc[1:-1].to_numpy()   
    if sheet_name != "Fig. 4": 
        l_down = df["L-down"].iloc[1:-1].to_numpy()    
        l_up = df["L-up"].iloc[1:-1].to_numpy()    
        l_star = l_up + l_down
    else:
        l_star = df["L-*"].iloc[1:-1].to_numpy()  
    q = df["Q"].iloc[1:-1].to_numpy()   

    K_down_values, K_up_values, L_star_values, R_N_values, time_shift = radiation_model(scenario[sheet_name], time)
    
    plt.figure(figsize=(10, 6))
    
    plt.plot(time, K_down_values, label='K-down', color='red')
    plt.plot(time, K_up_values, label='K-up', color='cyan')
    plt.plot(time, L_star_values, label='L-up+L-down', color='green')
    plt.plot(time, R_N_values, label='Q', color='purple')
    plt.plot(t, k_down, label="K-down (experimental)", linestyle="--", color='red')
    plt.plot(t, k_up, label="K-up (experimental)", linestyle="--", color='cyan')
    plt.plot(t, l_star, label="L-up + L-down (experimental)", linestyle="--", color='green')
    plt.plot(t, q, label="Q (experimental)", linestyle="--", color='purple')
    
    plt.title(f"Radiative Fluxes - {scenario[sheet_name]['title']}")
    plt.xlabel("Time [hrs]")
    plt.ylabel("Radiative Flux [W/m²]")
    plt.xlim(0,24)
    plt.legend()
    plt.grid()
    plt.tight_layout()
    
    plt.savefig(f"plot_{sheet_name}.png") 
    plt.show() 
 
    
 
    
 













    
 
#================================================================================    
#================================================================================ 
#================================================================================
#================================================================================    
#============================= Uncertainties ====================================




















#---------------- Definition of differential fuctions for uncertainties ---------------

def sinus_psi(phi, delta_s, t, lambda_e):
    """
    param phi: latitude of the location
    param delta_s: solar declination angle
    param t: time of the day in hours
    param lambda_e: longitude of the location
    """
    sin_psi = math.sin(phi) * math.sin(delta_s) - math.cos(phi) * math.cos(delta_s) * math.cos((math.pi * t / 12) + lambda_e)
    return sin_psi

def Delta_solar_declination_angle(data, hour=0, minute=0, second=0):
    """
    Calculate the change on solar declination angle.
    """
    year = data["year"]
    month = data["month"]
    day = data["day"]
    
    phi_r = 0.409 
    jr = 173  
    jt = time_of_year_from_date(datetime(year, month, day, hour, minute, second)) 
    J = 365
    
    D_J = data["D_J"]
    D_jt = data["D_j_t"]
    D_jr = data["D_j_r"]
    
    D_delta_s = abs(phi_r*2.0*np.pi/J*np.sin(2*np.pi*(jt-jr)/J)*D_jt) + abs(phi_r*2.0*np.pi/J*np.sin(2.0*np.pi*(jt-jr)/J)*D_jr) + abs(phi_r*2.0*np.pi*(jt-jr)/(J**2)*np.sin(2.0*np.pi*(jt-jr)/J)*D_J)
    
    return D_delta_s

def Delta_sinus_psi(t, data):
    """
    param phi: latitude of the location
    param delta_s: solar declination angle
    param t: time of the day in hours
    param lambda_e: longitude of the location
    """
    phi = data["phi"]
    lambda_e = data["lambda_e"]
    delta_s = solar_declination_angle(data["year"], data["month"], data["day"])
    D_delta_s = Delta_solar_declination_angle(data)
    
    Dt = data["Dt"]
    
    D_sin_psi = abs(math.cos(phi)*math.cos(math.pi*t/12+lambda_e)*math.sin(delta_s)*D_delta_s) + abs(math.cos(phi)*math.cos(delta_s)*math.pi/12*math.sin(math.pi*t/12+lambda_e)*Dt)
    return D_sin_psi


def Delta_Tau_K(t, data):
    phi = data["phi"]
    lambda_e = data["lambda_e"]
    delta_s = solar_declination_angle(data["year"], data["month"], data["day"])  
    sin_psi = sinus_psi(phi, delta_s, t, lambda_e)
    
    D_sin_psi = Delta_sinus_psi(t, data)
    sigma_CH = data["sigma_C_H"](t)
    sigma_CL = data["sigma_C_L"](t)
    sigma_CM = data["sigma_C_M"](t)
    
    D_CH = data["D_sigma_C_H"]
    D_CL = data["D_sigma_C_L"]
    D_CM = data["D_sigma_C_M"]
    
    D_TK = abs(0.2*(1-0.4*sigma_CH)*(1-0.7*sigma_CM)*(1-0.4*sigma_CL)*D_sin_psi) + abs(0.4*(0.6+0.2*sin_psi)*(1-0.7*sigma_CM)*(1-0.4*sigma_CL)*D_CH) + abs(0.7*(0.6+0.2*sin_psi)*(1-0.4*sigma_CH)*(1-0.4*sigma_CL)*D_CM) + abs(0.4*(0.6+0.2*sin_psi)*(1-0.7*sigma_CM)*(1-0.4*sigma_CH)*D_CL)
    return D_TK

def Delta_K_down(t, data):
    phi = data["phi"]
    lambda_e = data["lambda_e"]
    delta_s = solar_declination_angle(data["year"], data["month"], data["day"])  
    sin_psi = sinus_psi(phi, delta_s, t, lambda_e)
    sigma_CH = data["sigma_C_H"](t)
    sigma_CL = data["sigma_C_L"](t)
    sigma_CM = data["sigma_C_M"](t)
    
    TK = (0.6 + 0.2 * sin_psi) * (1 - 0.4 * sigma_CH) * (1 - 0.7 * sigma_CM) * (1 - 0.4 * sigma_CL)
    D_TK = Delta_Tau_K(t, data)
    D_sin_psi = Delta_sinus_psi(t, data)
    D_S0 = data["D_S0"]
    S0 = 1370 # W/m^2 (solar constant)
    
    D_Kd = abs(TK*sin_psi*D_S0) + abs(S0*sin_psi*D_TK) + abs(S0*TK*D_sin_psi)
    return D_Kd

def Delta_K_up(t, data):
    sigma_CH = data["sigma_C_H"](t)
    sigma_CL = data["sigma_C_L"](t)
    sigma_CM = data["sigma_C_M"](t)
    delta_s = solar_declination_angle(data["year"], data["month"], data["day"])  
    psi = solar_elevation_angle(data["phi"], delta_s, t, data["lambda_e"])
    
    alpha = data["albedo"]
    D_alpha = data["D_alpha"]
    
    Kd = incoming_shortwave_radiation(psi, sigma_CH, sigma_CM, sigma_CL)
    D_Kd = Delta_K_down(t, data)
    return abs(Kd*D_alpha) + abs(alpha*D_Kd)

def Delta_L(data):
    D_CH = data["D_sigma_C_H"]
    D_CL = data["D_sigma_C_L"]
    D_CM = data["D_sigma_C_M"]
    D_L = abs(97.28*0.1*D_CH) + abs(97.28*0.3*D_CM) + abs(97.28*0.6*D_CL)
    return D_L

def Delta_Net(t, data):
    D_Kd = Delta_K_down(t, data)
    D_Ku = Delta_K_up(t, data)
    D_L = Delta_L(data)
    return D_Kd + D_Ku + D_L

def Delta_radiation_model(data,time):
    Delta_K_down_values,  Delta_K_up_values, Delta_L_star_values,  Delta_R_N_values = [], [], [], []
    utc = UTC(data["lambda_e"])
    time_shift = [t_shift(t, utc) for t in time]
    # plt.figure(figsize=(10, 6))
    # plt.plot(time, sigma_C_H, label='sigma_C_H')
    # plt.plot(time, sigma_C_M, label='sigma_C_M')
    # plt.plot(time, sigma_C_L, label='sigma_C_L')
    for i, t in enumerate(time_shift):

        D_K_down = Delta_K_down(t, data)
        D_K_up =Delta_K_up(t, data)
        D_L_star = Delta_L(data)
        D_R_N = Delta_Net(t, data)

        Delta_K_down_values.append(D_K_down)
        Delta_K_up_values.append(D_K_up)
        Delta_L_star_values.append(D_L_star)
        Delta_R_N_values.append(D_R_N)

    return Delta_K_down_values, Delta_K_up_values, Delta_L_star_values, Delta_R_N_values, time_shift




#=========================== Main program part ========================




for sheet_name in excel_data.sheet_names:
    df = pd.read_excel(file_path, sheet_name=sheet_name)
    
    t = df["t "].iloc[1:-1].to_numpy()    
    k_down = df["K-down"].iloc[1:-1].to_numpy()  
    k_up = df["K-up"].iloc[1:-1].to_numpy()   
    if sheet_name != "Fig. 4": 
        l_down = df["L-down"].iloc[1:-1].to_numpy()    
        l_up = df["L-up"].iloc[1:-1].to_numpy()    
        l_star = l_up + l_down
    else:
        l_star = df["L-*"].iloc[1:-1].to_numpy()  
    q = df["Q"].iloc[1:-1].to_numpy()   

    K_down_values, K_up_values, L_star_values, R_N_values, time_shift = radiation_model(scenario[sheet_name], time)
    
    Delta_K_down_values, Delta_K_up_values, Delta_L_star_values, Delta_R_N_values, time_shift = Delta_radiation_model(scenario[sheet_name], time)
    
    K_d_sup = []
    K_d_inf = []
    K_u_sup = []
    K_u_inf = []
    L_sup = []
    L_inf = []
    Net_sup = []
    Net_inf = []
    
    
    for i in range(len(K_down_values)):
        K_d_sup.append(K_down_values[i] + 0.5*Delta_K_down_values[i])
        K_d_inf.append(K_down_values[i] - 0.5*Delta_K_down_values[i])
        K_u_sup.append(K_up_values[i] + 0.5*Delta_K_up_values[i])
        K_u_inf.append(K_up_values[i] - 0.5*Delta_K_up_values[i])
        L_sup.append(L_star_values[i] + 0.5*Delta_L_star_values[i])
        L_inf.append(L_star_values[i] - 0.5*Delta_L_star_values[i])
        Net_sup.append(R_N_values[i] + 0.5*Delta_R_N_values[i])
        Net_inf.append(R_N_values[i] - 0.5*Delta_R_N_values[i])
        
    plt.figure(figsize=(10, 6))
        
    plt.plot(time, K_down_values, label='K-down', color='red')
    plt.plot(time, K_up_values, label='K-up', color='cyan')
    plt.plot(time, L_star_values, label='L-up+L-down', color='green')
    plt.plot(time, R_N_values, label='Q', color='purple')
    
    plt.plot(time, K_d_sup, color='red', linestyle="--")
    plt.plot(time, K_u_sup, color='cyan', linestyle="--")
    plt.plot(time, L_sup, color='green', linestyle="--")
    plt.plot(time, Net_sup, color='purple', linestyle="--")
    
    plt.plot(time, K_d_inf, color='red', linestyle="--")
    plt.plot(time, K_u_inf, color='cyan', linestyle="--")
    plt.plot(time, L_inf, color='green', linestyle="--")
    plt.plot(time, Net_inf, color='purple', linestyle="--")
    

    plt.title(f"Radiative Fluxes - {scenario[sheet_name]['title']}")
    plt.xlabel("Time [hrs]")
    plt.ylabel("Radiative Flux [W/m²]")
    plt.xlim(0,24)
    plt.legend()
    plt.grid()
    plt.tight_layout()
    
    plt.savefig(f"uncertainties_{sheet_name}.png") 
    plt.show() 
    
 
    