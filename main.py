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


time = np.arange(0, 24, 0.05)  # Time in hours

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