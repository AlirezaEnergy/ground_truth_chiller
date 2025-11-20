from CoolProp.CoolProp import PropsSI
from CoolProp.HumidAirProp import HAPropsSI
import math



### supply fan model
def fan(supply_air_mass_flowrate, enthalpy_mix, humidity_ratio_mix):
    
    # ---- constants
    fan_rated_power = 560          ## [W]
    fan_rated_mass_flowrate = 1.33 ## [kg/s]
    barometric_pressure = 101325.0 # [Pa] equal to 1 atm
    c1 = 0.0015302446              # c1-c5: from energyplus variable speed fan
    c2 = 0.0052080574
    c3 = 1.1086242
    c4 = -0.11635563
    c5 = 0.000
    
    # ---- fan inlet conditions
    enthalpy_fan_in = enthalpy_mix             # [C] fan inlet temperature
    humidity_ratio_in_fan = humidity_ratio_mix # [kgWater/kgDryAir] fan inlet humidity ratio
    
    # ---- calculate fan power consumption
    FF = supply_air_mass_flowrate / fan_rated_mass_flowrate # [0-1] fan flow fraction
    PLF = c1 + c2*FF + c3*FF**2 + c4*FF**3 + c5*FF**4       # [0-1] fan part load factor
    fan_power = fan_rated_power * PLF                       # [kW] fan power consumption
    
    # ---- calculate fan outlet condition
    """
    since fan motor is located inside the air stream based on E+ and ASHRAE fundamentals 2021 page 509
    all power consumed by the fan is dissipated to the air stream
    """
    enthalpy_out_fan = enthalpy_fan_in + fan_power/supply_air_mass_flowrate                                                   # [J/kg]
    humidity_ratio_out_fan = humidity_ratio_in_fan                                                                            # [kgWater/kgDryAir]
    Tfan_out = HAPropsSI('T', 'H', enthalpy_out_fan, 'W', humidity_ratio_out_fan, 'P', barometric_pressure)-273.15            # [C]
    Tfan_out_dewpoint = HAPropsSI('Tdp', 'H', enthalpy_out_fan, 'W', humidity_ratio_out_fan, 'P', barometric_pressure)-273.15 # [C]
    relative_humidity_out_fan = HAPropsSI('RH', 'T', Tfan_out+273.15, 'W', humidity_ratio_out_fan, 'P', barometric_pressure)  # [0-1]
    
    return Tfan_out, humidity_ratio_out_fan, relative_humidity_out_fan, Tfan_out_dewpoint, enthalpy_out_fan, fan_power



def cooling_coil(supply_air_mass_flowrate, 
                 supply_air_temp, 
                 supply_air_saturation_humidity_ratio, 
                 Tfan_out, 
                 Tfan_out_dewpoint, 
                 humidity_ratio_out_fan, 
                 enthalpy_out_fan):
    
    # ---- constants
    barometric_pressure = 101325.0 # [Pa] equal to 1 atm
    UA = 3                         # [kW/C] cooling coil UA value
    
    # ---- coil inlet condition
    Tcoil_in               = Tfan_out               # [C] cooling coil inlet temperature
    Tcoil_in_dewpoint      = Tfan_out_dewpoint      # [C] cooling coil inlet dewpoint temperature
    humidity_ratio_in_coil = humidity_ratio_out_fan # [kgWater/kgDryAir] cooling coil inlet humidity ratio
    enthalpy_in_coil       = enthalpy_out_fan       # [J/kgDryAir]
    
    # ---- coil outlet condition
    Tcoil_out = supply_air_temp  # [C] cooling coil outlet initial temperature
    
    c_p_air = HAPropsSI('cp', 'T', Tcoil_in+273.15, 'W', humidity_ratio_in_coil, 'P', barometric_pressure)/1000 # [kJ/kgC] coil inlet specific heat
    BF = math.exp( -UA/(supply_air_mass_flowrate*c_p_air) ) # [0-1] cooling coil bypass factor
    
    # ---- calculating coil outlet humidity ratio
    
    Ts = (Tcoil_out - BF*Tcoil_in)/(1-BF) # [C] coil surface temperature

    # ---- finding coil outlet humidity ratio
    if Ts < Tcoil_in_dewpoint: # condensation occurs

        humidity_ratio_surface_coil = HAPropsSI('W', 'T', Ts+273.15, 'RH', 1.0, 'P',barometric_pressure)#[kgWater/kgDryAir]cooling coil surface humidityratio
        humidity_ratio_out_coil = BF*humidity_ratio_in_coil +(1-BF)*humidity_ratio_surface_coil         #[kgWater/kgDryAir]cooling coil outlet humidity ratio
        
    else: # Ts >= Tcoil_in_dewpoint and condensation does not occur

        humidity_ratio_out_coil = humidity_ratio_in_coil  # [kgWater/kgDryAir]
    
    counter = 0
    while humidity_ratio_out_coil > supply_air_saturation_humidity_ratio and counter < 200:
        
        Tcoil_out += 0.1
        supply_air_saturation_humidity_ratio = HAPropsSI('W', 'T', Tcoil_out+273.15, 'RH', 1.0, 'P',barometric_pressure) # [kgWater/kgDryAir]
        
        Ts = (Tcoil_out - BF*Tcoil_in)/(1-BF) # [C] coil surface temperature

        # ---- finding coil outlet humidity ratio

        if Ts < Tcoil_in_dewpoint: # condensation occurs

            humidity_ratio_surface_coil = HAPropsSI('W', 'T', Ts+273.15, 'RH', 1.0, 'P',barometric_pressure)#[kgWater/kgDryAir]cooling coil surface humidityratio
            humidity_ratio_out_coil = BF*humidity_ratio_in_coil +(1-BF)*humidity_ratio_surface_coil         #[kgWater/kgDryAir]cooling coil outlet humidity ratio
            
        else: # Ts >= Tcoil_in_dewpoint and condensation does not occur

            humidity_ratio_out_coil = humidity_ratio_in_coil  # [kgWater/kgDryAir]
            
        counter += 1
    
    # ---- calculating mass flowrate and temperature of the condensate
    """
    Temperature of the condensate leaving the system is assumed to be equal to the wet-bulb temperature at coil outlet condition.
    This assumption is based on ASHRAE Principles of Heating Ventilation and Air Conditioning (2013) page 481
    """
    m_dot_condensate = supply_air_mass_flowrate*(humidity_ratio_in_coil - humidity_ratio_out_coil)                 # [kgWater/s]
    T_condensate = HAPropsSI('Twb', 'T', Tcoil_out+273.15, 'W', humidity_ratio_out_coil, 'P', barometric_pressure) # [K]
    
    # ---- calculating condensate water enthalpy and energy
    """
    It is assumed that leaving condensate water is saturated (i.e., its quality is 0).
    This assumption is based on ASHRAE Principles of Heating Ventilation and Air Conditioning (2013) page 481
    """
    enthalpy_condensate = PropsSI("H", "T", T_condensate, "Q", 0, "water") # [J/kgWater] 273.15 is not added since T_condensate is in K
    Q_condensate = m_dot_condensate*enthalpy_condensate/1000 # [kW] divided by 1000 to turn W to kW
    
    # ---- calculating coil outlet enthalpy and relative humidity
    enthalpy_out_coil = HAPropsSI('H', 'T', Tcoil_out+273.15, 'W', humidity_ratio_out_coil, 'P', barometric_pressure)           # [J/kgDryAir]
    relative_humidity_out_coil = HAPropsSI('RH', 'T', Tcoil_out+273.15, 'W', humidity_ratio_out_coil, 'P', barometric_pressure) # [0-1]
    
    # ---- performing first law analysis (i.e., energy balance) on the entire cooling coil
    """
    total_cooling_load = sensible_load + latent_load - condensate_energy
    sensible_load + latent_load = delta_enthalpy of inlet and outlet
    """
    Q_sensible_and_latent = supply_air_mass_flowrate*(enthalpy_in_coil - enthalpy_out_coil)/1000 # [kW] total energy absorbed from the moist air (/1000:W=>kW)
    total_cooling_load = Q_sensible_and_latent - Q_condensate                                    # [kW] total energy transfered to the refrigerant
    
    return total_cooling_load, humidity_ratio_out_coil, relative_humidity_out_coil, Tcoil_out, Ts


def chiller_on_off_status(supply_air_mass_flowrate, 
                          supply_air_temp, 
                          supply_air_saturation_humidity_ratio,
                          enthalpy_mix, 
                          humidity_ratio_mix):
    
    barometric_pressure = 101325.0   # [Pa] equal to 1 atm
    chiller_min_acceptable_load = 0  ## [kW]
    chiller_max_acceptable_load = 25 ## [kW]
    
    chiller_on = True
    
    Tfan_out, humidity_ratio_out_fan, relative_humidity_out_fan, Tfan_out_dewpoint, enthalpy_out_fan, fan_power = fan(supply_air_mass_flowrate, 
                                                                                                                      enthalpy_mix, 
                                                                                                                      humidity_ratio_mix)
    
    Qe, supply_air_humidity_ratio, supply_air_relative_humidity, supply_air_temp, T_surface_coil = cooling_coil(supply_air_mass_flowrate, 
                                                                                                                supply_air_temp, 
                                                                                                                supply_air_saturation_humidity_ratio, 
                                                                                                                Tfan_out, 
                                                                                                                Tfan_out_dewpoint, 
                                                                                                                humidity_ratio_out_fan, 
                                                                                                                enthalpy_out_fan)
    
    if Qe < chiller_min_acceptable_load: # chiller wont turn on
        
        chiller_on = False
    
    elif Qe > chiller_max_acceptable_load: # chiller operates at its max capacity
        
        counter = 0
        while Qe > chiller_max_acceptable_load and counter < 200:
    
            supply_air_temp += 0.1 # [C], increase supply air temp to reduce the cooling load on coil
            supply_air_saturation_humidity_ratio = HAPropsSI('W', 'T', supply_air_temp+273.15, 'RH', 1.0, 'P',barometric_pressure) # [kgWater/kgDryAir]
            
            Qe, supply_air_humidity_ratio, supply_air_relative_humidity, supply_air_temp, T_surface_coil = cooling_coil(supply_air_mass_flowrate, 
                                                                                                                        supply_air_temp, 
                                                                                                                        supply_air_saturation_humidity_ratio, 
                                                                                                                        Tfan_out, 
                                                                                                                        Tfan_out_dewpoint, 
                                                                                                                        humidity_ratio_out_fan, 
                                                                                                                        enthalpy_out_fan)
            counter += 1
    
    return chiller_on, Qe, supply_air_humidity_ratio, supply_air_relative_humidity, T_surface_coil, Tfan_out, humidity_ratio_out_fan, relative_humidity_out_fan, Tfan_out_dewpoint, enthalpy_out_fan, fan_power, supply_air_temp



def condenser_fan(outdoor_temp, Qc):
    
    fan_rated_power = 1000       ## [W]
    fan_rated_mass_flowrate = 2  ## [kg/s]
    c1 = 0.0015302446            # c1-c5: from energyplus variable speed fan
    c2 = 0.0052080574
    c3 = 1.1086242
    c4 = -0.11635563
    c5 = 0.000
    DeltaT = 8  # air temperature uplift across the condenser
    
    cooling_air_mass_flow_rate = min(Qc/(1.005*DeltaT), fan_rated_mass_flowrate)
    T_air_out = outdoor_temp + Qc/(cooling_air_mass_flow_rate*1.005)
        
    FF = cooling_air_mass_flow_rate / fan_rated_mass_flowrate # [0-1] fan flow fraction
    PLF = c1 + c2*FF + c3*FF**2 + c4*FF**3 + c5*FF**4         # [0-1] fan part load factor
    fan_power = fan_rated_power * PLF                         # [kW] fan power consumption
    
    return cooling_air_mass_flow_rate, fan_power, T_air_out



def ground_truth_chiller(outdoor_temp, 
                         supply_air_mass_flowrate,
                         supply_air_temp,
                         supply_air_saturation_humidity_ratio,
                         Tmix,
                         Tdew_mix,
                         humidity_ratio_mix,
                         enthalpy_mix):
    
    # ---- constants
    
    T_superheat_evaporator = 5   # [C]     evaporator, temperature increase in the superheater section of the evaporator
    T_subcooling_condenser = 2   # [C]     temperature decrease after the two-phase region of the condenser
    A = 0                        # [C/C^2] coef. of T_superheat_evaporator^2 in E+ quadratic curve for calculating evaporating temperature
    B = 0.843                    # [-]     coef. of T_superheat_evaporator in E+ quadratic curve for calculating evaporating temperature
    C = 0                        # [C]     constant coef. of E+ quadratic curve for calculating evaporating temperature
    
    DeltaT_cooling_air = 12      # [C] condenser,  condensing temp = air temp + DeltaT_cooling_air
    Tmin_cond = 21               # [C] condenser minimum temperature
    min_eva_cond_T_diff = 10     # [C] minimum temperature difference between the evaporator and condenser
    Pressure_drop_superheater = 5/100 # 5% pressure drop in the superheat region
    eta_mechanical = 0.90        # [0-1] compressor mechanical efficiency
    W_max_motor = 7.15           ## [kW] max shaft power
    eta_max_motor = 0.95         # [0-1] max motor electrical efficiency
    eta_motor_slope = 10         # [-]
        
    Vd = 192.2e-6                # [m3/revolution]
    Av = 0.00683154              # coefficients of volumetric efficiency
    Bv = -0.06876312
    Cv = 1.00797235
    
    # divide the supply_air_mass_flowrate equally between the units
    supply_air_mass_flowrate = supply_air_mass_flowrate/4
    
    chiller_on, Qe, supply_air_humidity_ratio, supply_air_relative_humidity, T_surface_coil, Tfan_out, humidity_ratio_out_fan, relative_humidity_out_fan, Tfan_out_dewpoint, enthalpy_out_fan, fan_power, supply_air_temp = chiller_on_off_status(supply_air_mass_flowrate, 
                                                                                                                                                                                                                                                  supply_air_temp, 
                                                                                                                                                                                                                                                  supply_air_saturation_humidity_ratio,
                                                                                                                                                                                                                                                  enthalpy_mix, 
                                                                                                                                                                                                                                                  humidity_ratio_mix)
    
    if chiller_on:
        
        # ---- find evaporating temperature
        
        Te = T_surface_coil - (A*T_superheat_evaporator**2 + B*T_superheat_evaporator + C)  # [C] E+ equation
        
        ##### evaporator
        
        T1 = Te                                                  # [C]
        P1 = PropsSI("P", "T", T1+273.15, "Q", 1, "R134a")/1000  # [kPa]
        
        ##### superheater
        
        T1_2 = T1 + T_superheat_evaporator                                   # [C]
        P1_2 = P1*(1-Pressure_drop_superheater)                              # [kPa]
        h1_2 = PropsSI("H", "T", T1_2+273.15, "P", P1_2*1000, "R134a")/1000  # [kJ/kg]
        s1_2 = PropsSI("S", "T", T1_2+273.15, "P", P1_2*1000, "R134a")       # [J/kgC]
        rho1_2 = PropsSI("D", "H", h1_2*1000, "P", P1_2*1000, "R134a")       # [kg/m3]
        
        ##### condenser
        
        Tmin_cond = max(Tmin_cond, Te + min_eva_cond_T_diff)
        Tc = max(outdoor_temp + DeltaT_cooling_air, Tmin_cond)     # [C] condensing temperature
        T3 = Tc                                                    # [C]
        P3 = PropsSI("P", "T", T3+273.15, "Q", 0, "R134a")/1000    # [kPa]
        
        ##### subcooler
        
        P3_2 = P3
        T3_2 = T3 - T_subcooling_condenser
        h3_2 = PropsSI("H", "T", T3_2+273.15, "P", P3_2*1000, "R134a")/1000  # [kJ/kg]
        rho3_2 = PropsSI("D", "H", h3_2*1000, "P", P3_2*1000, "R134a")       # [kg/m3] coolprop prefers using H instead of T here (T gives an error)
        
        ##### expansion valve
        
        P4 = P1                                                   # [kPa]
        h4 = h3_2                                                 # [kJ/kg]
        rho4 = PropsSI("D", "H", h4*1000, "P", P4*1000, "R134a")  # [kg/m3]
        
        mDotRef = Qe / (h1_2 - h4)  # [kg/s] refrigerant mass flow rate
        
        pressure_difference = (P3_2 - P4)*1000                         # [Pa] pressure drop in the electronic expansion valve (EEV)
        flow_coefficient = 0.02005 * (rho3_2)**0.5 + 0.634 * (1/rho4)  # [-]
        
        A_EEV = mDotRef / (flow_coefficient * abs(rho3_2 * pressure_difference)**0.5)  # [m2] EEV flow area
        
        ##### compressor
        
        P2 = P3       # [kPa]
        Pr = P2/P1_2  # [-]
        
        eff_s = 0.11229 + 0.34692*Pr - 0.068534*Pr**2 + 0.0046362*Pr**3  # [0-1] from paper 96
        
        s2s = s1_2                                                       # [J/kgC]
        h2s = PropsSI("H", "P", P2*1000, "S", s2s, "R134a")/1000         # [kJ/kg]
        h2 = h1_2 + (h2s-h1_2)/eff_s                                     # [kJ/kg]
        
        eff_v = Av*Pr**2 + Bv*Pr + Cv                # [0-1]  volumetric efficiency
        V_dot_actual = mDotRef/rho1_2                # [m3/s] actual volumetric flowrate
        V_dot_theory = V_dot_actual/eff_v            # [m3/s] theoritical volumetric flowrate
        compressor_frequency = (V_dot_theory*60)/Vd  # [rpm]  compressor frequency or revolution speed
        
        W_cycle = mDotRef*(h2 - h1_2)     # [kW]
        W_shaft = W_cycle/eta_mechanical  # [kW]
        PLR_motor = W_shaft / W_max_motor # [0-1]
        
        eta_motor = eta_max_motor * ( 2/(1+math.exp(-eta_motor_slope*PLR_motor)) - 1 ) # [0-1]
        W_motor = W_shaft / eta_motor                                                  # [kW]
        
        
        Qc = mDotRef*(h2-h3_2)  # [kW]
        
        cooling_air_mass_flow_rate, cooling_fan_power, cooling_T_air_out = condenser_fan(outdoor_temp, Qc)
        
        first_law_error = abs(W_cycle + Qe - Qc)  # [kW]
        

        output_dict = {
                       "Qe": Qe*4,
                       "Te": Te,
                       "Tc": Tc,
                       "mDotRef": mDotRef*4,
                       "Pr": Pr,
                       "eff_s": eff_s,
                       "W_motor": W_motor*4,
                       "first_law_error": first_law_error,
                       "supply_air_humidity_ratio": supply_air_humidity_ratio,
                       "supply_air_relative_humidity": supply_air_relative_humidity,
                       "Tfan_out": Tfan_out,
                       "fan_power": fan_power*4,
                       "A_EEV": A_EEV,
                       "PLR_motor": PLR_motor,
                       "compressor_frequency":compressor_frequency,
                       "eff_v":eff_v,
                       "cooling_air_mass_flow_rate": cooling_air_mass_flow_rate*4,
                       "cooling_T_air_out": cooling_T_air_out,
                       "cooling_fan_power": cooling_fan_power*4,
                       "supply_air_temp": supply_air_temp
                      }
        
    else:
        
        output_dict = {
                       "Qe": -10,
                       "Te": -10,
                       "Tc": -10,
                       "mDotRef": -10,
                       "Pr": -10,
                       "eff_s": -10,
                       "W_motor": -10,
                       "first_law_error": -10,
                       "supply_air_humidity_ratio": humidity_ratio_out_fan,
                       "supply_air_relative_humidity": relative_humidity_out_fan,
                       "Tfan_out": Tfan_out,
                       "fan_power": fan_power*4,
                       "A_EEV": -10,
                       "PLR_motor": -10,
                       "compressor_frequency":-10,
                       "eff_v":-10,
                       "cooling_air_mass_flow_rate": -10,
                       "cooling_T_air_out": -10,
                       "cooling_fan_power": -10,
                       "supply_air_temp": supply_air_temp
                      }
    
    return output_dict









