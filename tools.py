#---------------------------------------------------- FUNCOES TERMODINAMICAS ------------------------------------------#
import numpy as np
import scipy.optimize as sp_opt

#---------------------------------------------- CALCULO DE PRESSAO DE SATURACAO ---------------------------------------#
def Psat(T_K: float, C1, C2, C3, C4, C5):
    ln_Psat = C1 + (C2 / T_K) + C3 * np.log(T_K) + C4 * (T_K ** C5)
    return np.exp(ln_Psat)

#-------------------------------------------- CALCULO DE COEFICIENTES DE ATIVIDADE -------------------------------------#
def Coef_atv(modelo, ctes_thermo, T_K):
    x1 = np.linspace(0.00, 1.00, 101)        
    x2 = 1.0 - x1

    match modelo:
        case 'Ideal':
            return np.ones_like(x1), np.ones_like(x1)
        
        case 'NRTL':
            a12 = ctes_thermo[modelo].get('a12', 0.0)
            a21 = ctes_thermo[modelo].get('a21', 0.0)
            b12 = ctes_thermo[modelo].get('b12', 0.0)
            b21 = ctes_thermo[modelo].get('b21', 0.0)
            c12 = ctes_thermo[modelo].get('c12', 0.30)
            c21 = ctes_thermo[modelo].get('c21', 0.30)

            tau12 = a12 + (b12 / T_K)
            tau21 = a21 + (b21 / T_K)

            G12 = np.exp(-c12 * tau12)
            G21 = np.exp(-c21 * tau21)

            term1_gamma1 = tau21 * ((G21 / (x1 + G21 * x2)) ** 2)
            term2_gamma1 = (tau12 * G12) / ((G12 * x1 + x2) ** 2)
            ln_gamma1 = (x2 ** 2) * (term1_gamma1 + term2_gamma1)

            term1_gamma2 = tau12 * ((G12 / (G12 * x1 + x2)) ** 2)
            term2_gamma2 = (tau21 * G21) / ((x1 + G21 * x2) ** 2)
            ln_gamma2 = (x1 ** 2) * (term1_gamma2 + term2_gamma2)

            return np.exp(ln_gamma1), np.exp(ln_gamma2)

        case 'Wilson':
            a12 = ctes_thermo[modelo].get('a12', 0.0)
            a21 = ctes_thermo[modelo].get('a21', 0.0)
            b12 = ctes_thermo[modelo].get('b12', 0.0)
            b21 = ctes_thermo[modelo].get('b21', 0.0)
            c12 = ctes_thermo[modelo].get('c12', 0.0)
            c21 = ctes_thermo[modelo].get('c21', 0.0)

            ln_A12 = a12 + (b12 / T_K) + (c12 * np.log(T_K))
            ln_A21 = a21 + (b21 / T_K) + (c21 * np.log(T_K))

            L12 = np.exp(ln_A12)
            L21 = np.exp(ln_A21)

            denom1 = x1 + x2 * L12
            denom2 = x1 * L21 + x2

            ln_gamma1 = 1.0 - np.log(denom1) - (x1 / denom1) - (x2 * L21 / denom2)
            ln_gamma2 = 1.0 - np.log(denom2) - (x2 / denom2) - (x1 * L12 / denom1)

            return np.exp(ln_gamma1), np.exp(ln_gamma2)
        
        case _:
            raise ValueError(f"Modelo {modelo} não suportado.")


#-------------------------------------------- CALCULO DE TEMPERATURA DE SATURACAO --------------------------------------#
def temperatura_sat_componente(P_total_Pa, comp_dict):
    
    def erro_sat(T):
        return Psat(T, **comp_dict['Antoine']) - P_total_Pa
    
    Tmin = comp_dict['Limites']['Tmin']
    Tmax = comp_dict['Limites']['Tmax']
    try:
        return float(sp_opt.brentq(erro_sat, Tmin, Tmax))
    except ValueError:
        P_tmin = Psat(Tmin, **comp_dict['Antoine'])
        return float(Tmin if P_total_Pa < P_tmin else Tmax)

#---------------------------------------------- CALCULO DE TEMPERATURA DE BOLHA ----------------------------------------#
def temperatura_bolha(x1_val, P_total_Pa, modelo, ctes_thermo, db_compostos, nome_v, nome_p):
  
    x1 = float(np.clip(x1_val, 1e-6, 1.0 - 1e-6))
    x2 = 1.0 - x1

    comp_v = db_compostos[nome_v]
    comp_p = db_compostos[nome_p]

    T_sat_v = temperatura_sat_componente(P_total_Pa, comp_v)
    T_sat_p = temperatura_sat_componente(P_total_Pa, comp_p)

    T_low = max(150.0, min(comp_v['Limites']['Tmin'], comp_p['Limites']['Tmin']))
    T_high = min(650.0, max(comp_v['Limites']['Tmax'], comp_p['Limites']['Tmax']))

    def erro_pressao(T_k):
        Psat1 = Psat(T_k, **comp_v['Antoine'])
        Psat2 = Psat(T_k, **comp_p['Antoine'])

        g_v_arr, g_p_arr = Coef_atv(modelo, ctes_thermo, T_k)
        x_grid = np.linspace(0.0, 1.0, 101)
        gamma1 = float(np.interp(x1, x_grid, g_v_arr))
        gamma2 = float(np.interp(x1, x_grid, g_p_arr))

        P_calc = x1 * gamma1 * Psat1 + x2 * gamma2 * Psat2
        return P_calc - P_total_Pa

    grid_T = np.linspace(T_low, T_high, 120)
    erros = np.array([erro_pressao(t) for t in grid_T])
    
    sign_changes = np.where(np.diff(np.sign(erros)))[0]

    if len(sign_changes) > 0:
        idx = sign_changes[0]
        T_a, T_b = grid_T[idx], grid_T[idx + 1]
        try:
            T_sol = float(sp_opt.brentq(erro_pressao, T_a, T_b))
        except ValueError:
            T_sol = float(grid_T[np.argmin(np.abs(erros))])
    else:
        res = sp_opt.minimize_scalar(lambda t: abs(erro_pressao(t)), bounds=(T_low, T_high), method='bounded')
        T_sol = float(res.x)

    Psat1_sol = Psat(T_sol, **comp_v['Antoine'])
    g_v_arr, _ = Coef_atv(modelo, ctes_thermo, T_sol)
    gamma1_sol = float(np.interp(x1, np.linspace(0.0, 1.0, 101), g_v_arr))
    
    P1_parcial = x1 * gamma1_sol * Psat1_sol
    y1 = P1_parcial / P_total_Pa

    return float(T_sol), float(np.clip(y1, 0.0, 1.0))

#------------------------------------------ CALCULO DE ENTALPIA DE VAPORIZACAO -----------------------------------------#
def hvap_componente_watson(T_K, comp_dict):
    
    watson_data = comp_dict['Watson']
    Hvap_0 = watson_data['Hvap_0']
    Teb_0 = watson_data['Teb_0']
    Tc = watson_data['Tc']

    T_eval = min(T_K, Tc - 0.5)
    
    tr_atual = 1.0 - (T_eval / Tc)
    tr_ref = 1.0 - (Teb_0 / Tc)

    return Hvap_0 * ((tr_atual / tr_ref) ** 0.38)

def hvap_mistura_watson(T_K, x1, db_compostos, nome_v, nome_p):
    
    comp_v = db_compostos[nome_v]
    comp_p = db_compostos[nome_p]

    hvap1 = hvap_componente_watson(T_K, comp_v)
    hvap2 = hvap_componente_watson(T_K, comp_p)

    return x1 * hvap1 + (1.0 - x1) * hvap2
