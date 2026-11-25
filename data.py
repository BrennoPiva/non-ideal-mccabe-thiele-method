#---------------------------------------------------- DADOS E CONSTANTES -----------------------------------------------#
import streamlit as st

MODELOS = ['Ideal', 'NRTL', 'Wilson']

COMPOSTOS_PUROS = {
    'Etanol': {
        'Antoine': {'C1': 73.304, 'C2': -7122.3, 'C3': -7.1424, 'C4': 2.8853E-6, 'C5': 2},
        'Limites': {'Tmin': 159.05, 'P(Tmin)': 4.96E-4, 'Tmax': 514, 'P(Tmax)': 6.109E6},
        'Watson': {'Hvap_0': 38560.00, 'Teb_0': 351.5, 'Tc': 514.00}
    },
    'Água': {
        'Antoine': {'C1': 73.649, 'C2': -7258.2, 'C3': -7.3037, 'C4': 4.165E-6, 'C5': 2},
        'Limites': {'Tmin': 273.16, 'P(Tmin)': 6.11E2, 'Tmax': 425.56, 'P(Tmax)': 506625},
        'Watson': {'Hvap_0': 40660.00, 'Teb_0': 373.2, 'Tc': 647.1}
    },
    'Butanol': {
        'Antoine': {'C1': 106.295, 'C2': -9866.4, 'C3': -11.655, 'C4': 1.0832E-17, 'C5': 6},
        'Limites': {'Tmin': 183.85, 'P(Tmin)': 2.90E-4, 'Tmax': 563.1, 'P(Tmax)': 4.401E6},
        'Watson': {'Hvap_0': 43290.00, 'Teb_0': 390.6, 'Tc': 562.00}
    },
    'Acetona': {
        'Antoine': {'C1': 69.006, 'C2': -5599.6, 'C3': -7.0985, 'C4': 6.2237E-6, 'C5': 2},
        'Limites': {'Tmin': 178.45, 'P(Tmin)': 2.79, 'Tmax': 508.2, 'P(Tmax)': 4.709E6},
        'Watson': {'Hvap_0': 29100.00, 'Teb_0': 329.3, 'Tc': 508.00}
    },
    'Metanol': {
        'Antoine': {'C1': 82.718, 'C2': -6904.5, 'C3': -8.8622, 'C4': 7.4664E-6, 'C5': 2},
        'Limites': {'Tmin': 175.47, 'P(Tmin)': 1.11E-1, 'Tmax': 512.5, 'P(Tmax)': 8.146E6},
        'Watson': {'Hvap_0': 35210.00, 'Teb_0': 337.8, 'Tc': 513.00}
    }
}

MISTURAS = {
    'Metanol-Água': {
        'Volatil': 'Metanol',
        'Pesado': 'Água',
        'NRTL': {'a12': -0.693, 'a21': 2.7322, 'b12': 172.987, 'b21': -617.269, 'c12': 0.30, 'c21': 0.30},
        'Wilson': {'a12': -2.0302, 'a21': 0.049, 'b12': 448.679, 'b21': -21.9744}
    },
    'Acetona-Água': {
        'Volatil': 'Acetona',
        'Pesado': 'Água',
        'NRTL': {'a12': 6.3981, 'a21': 0.0544, 'b12': -1808.99, 'b21': 419.972, 'c12': 0.3, 'c21': 0.3},
        'Wilson': {'a12': -0.4374, 'a21': -9.7051, 'b12': -291.079, 'b21': 2719.21}
    },
    'Metanol-Etanol': {
        'Volatil': 'Metanol',
        'Pesado': 'Etanol',
        'NRTL': {'a12': 4.7119, 'a21': -2.3127, 'b12': -1162.29, 'b21': 483.844, 'c12': 0.3, 'c21': 0.3},
        'Wilson': {'a12': 1.5166, 'a21': -4.0246, 'b12': -288.173, 'b21': 991.718}
    },
    'Butanol-Água': {
        'Volatil': 'Água',
        'Pesado': 'Butanol',
        'NRTL': {'a12': 13.1102, 'a21': -2.0405, 'b12': -3338.95, 'b21': 763.869, 'c12': 0.3, 'c21': 0.3},
        'Wilson': {'a12': 1.1048, 'a21': 22.7777, 'b12': -588.733, 'b21': -10000}
    }
}

#--------------------------------------------- INICIALIZACAO DO BANCO EM SESSION --------------------------------------#
def inicializar_banco_dados():   
    if 'compostos_puros' not in st.session_state:
        st.session_state['compostos_puros'] = COMPOSTOS_PUROS.copy()
    if 'misturas' not in st.session_state:
        st.session_state['misturas'] = MISTURAS.copy()
    return st.session_state['compostos_puros'], st.session_state['misturas']