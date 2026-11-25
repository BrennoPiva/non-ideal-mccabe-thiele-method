#---------------------------------------------------- INTERFACE STREAMLIT ----------------------------------------------#
import streamlit as st
import tools as to
import backend as be
import data as db
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
import os


diretorio_atual = os.path.dirname(os.path.abspath(__file__))
caminho_unicamp = os.path.join(diretorio_atual, "logo_unicamp.jpg")
caminho_feq = os.path.join(diretorio_atual, "logo_feq.png")
#----------------------------------------------- CONFIGURACAO INICIAL DA PAGINA ----------------------------------------#
st.set_page_config(
    page_title="YuNno-Thermo",
    page_icon='♨️',
    layout='wide'
)

col_head1, col_head2 = st.columns([3, 1])

with col_head1:
    st.title("♨️ YuNno-Thermo")
    st.caption("Dimensionamento e simulação de destilação binária.")
    
    st.markdown(
        """
        <div style="padding: 10px 14px; border-radius: 6px; margin-bottom: 15px;">
            <p style="margin: 0; font-size: 0.85rem; color: #F5A623;">
                <b>Desenvolvido por:</b><br>
                <b>Brenno Fabricio Piva Ostroschi:</b> <i>@brennopiva | brenno.piva@gmail.com</i><br>
                <b>Yuri Stankunavicius:</b> <i>@yuri.stanks | yuri.ys.stank@gmail.com </i>
            </p>
        </div>
        """,
        unsafe_allow_html=True
    )

with col_head2:    
    col_l1, col_l2 = st.columns(2)
    with col_l1:
        if os.path.exists(caminho_unicamp):
            st.image(caminho_unicamp, use_container_width=True)
    with col_l2:
        if os.path.exists(caminho_feq):
            st.image(caminho_feq, use_container_width=True)

st.markdown("---")

#---------------------------------------------- BANCO DE DADOS DINAMICO DA SESSAO --------------------------------------#
COMPOSTOS_PUROS, MISTURAS = db.inicializar_banco_dados()

#---------------------------------------------- TRATAMENTO SEGURA DE ENTRADAS NUMERICAS --------------------------------#
def parse_num(val, default=0.0):
    if val is None:
        return float(default)
    if isinstance(val, (int, float)):
        return float(val)
    s = str(val).strip().replace(',', '.')
    if not s:
        return float(default)
    try:
        return float(s)
    except Exception:
        return float(default)


#---------------------------------------------- JANELA PARA ADICIONAR MISTURA -----------------------------------------#
@st.dialog("➕ Adicionar Nova Mistura Binária")
def modal_adicionar_mistura():
    import unicodedata

    def normalizar_texto(txt):
        if not txt:
            return ""
        txt_nfd = unicodedata.normalize('NFD', str(txt).strip().lower())
        return ''.join(c for c in txt_nfd if unicodedata.category(c) != 'Mn')

    
    mapa_compostos_norm = {
        normalizar_texto(nome): (nome, dados)
        for nome, dados in COMPOSTOS_PUROS.items()
    }

    def aplicar_autopreenchimento(prefixo_key, nome_digitado, flag_state_key):
        chave_norm = normalizar_texto(nome_digitado)
        ultimo_preenchido = st.session_state.get(flag_state_key, "")

        if chave_norm in mapa_compostos_norm and chave_norm != ultimo_preenchido:
            _, dados = mapa_compostos_norm[chave_norm]
            ant = dados['Antoine']
            lim = dados['Limites']
            wat = dados['Watson']

            st.session_state[f"c1_{prefixo_key}"] = str(ant['C1'])
            st.session_state[f"c2_{prefixo_key}"] = str(ant['C2'])
            st.session_state[f"c3_{prefixo_key}"] = str(ant['C3'])
            st.session_state[f"c4_{prefixo_key}"] = str(ant['C4'])
            st.session_state[f"c5_{prefixo_key}"] = str(ant['C5'])

            st.session_state[f"tmin_{prefixo_key}"] = str(lim['Tmin'])
            st.session_state[f"ptmin_{prefixo_key}"] = str(lim['P(Tmin)'])
            st.session_state[f"tmax_{prefixo_key}"] = str(lim['Tmax'])
            st.session_state[f"ptmax_{prefixo_key}"] = str(lim['P(Tmax)'])

            st.session_state[f"hvap_{prefixo_key}"] = str(wat['Hvap_0'])
            st.session_state[f"teb_{prefixo_key}"] = str(wat['Teb_0'])
            st.session_state[f"tc_{prefixo_key}"] = str(wat['Tc'])

            st.session_state[flag_state_key] = chave_norm
        elif chave_norm not in mapa_compostos_norm:
            st.session_state[flag_state_key] = ""

    st.caption(
        f"Insira os parâmetros da mistura. Compostos já cadastrados ({', '.join(COMPOSTOS_PUROS.keys())}) "
        "são preenchidos automaticamente ao digitar o nome."
    )

    tab1, tab2 = st.tabs(["Compostos Puros", "Interação Binária (NRTL / Wilson)"])

    with tab1:
        col_v, col_p = st.columns(2)
        with col_v:
            novo_v = st.text_input("Componente Volátil (Comp 1):", key="input_nome_v")
        with col_p:
            novo_p = st.text_input("Componente Pesado (Comp 2):", key="input_nome_p")

       
        aplicar_autopreenchimento("v", novo_v, "autopreenchido_comp1")
        aplicar_autopreenchimento("p", novo_p, "autopreenchido_comp2")

        if st.session_state.get("autopreenchido_comp1"):
            st.info(f"✅ Dados de **{novo_v.strip()}** carregados automaticamente do banco.")
        if st.session_state.get("autopreenchido_comp2"):
            st.info(f"✅ Dados de **{novo_p.strip()}** carregados automaticamente do banco.")

        st.markdown("---")
        st.markdown("**Antoine Estendida (DIPPR 101) (Comp 1):**")
        col_v_c1, col_v_c2, col_v_c3 = st.columns(3)
        c1_v = col_v_c1.text_input("C1", key="c1_v")
        c2_v = col_v_c2.text_input("C2", key="c2_v")
        c3_v = col_v_c3.text_input("C3", key="c3_v")

        col_v_c4, col_v_c5, col_v_tmin = st.columns(3)
        c4_v = col_v_c4.text_input("C4", key="c4_v")
        c5_v = col_v_c5.text_input("C5", key="c5_v")
        tmin_v = col_v_tmin.text_input("Tmin [K]", key="tmin_v")

        col_v_ptmin, col_v_tmax, col_v_ptmax = st.columns(3)
        ptmin_v = col_v_ptmin.text_input("P(Tmin) [Pa]", key="ptmin_v")
        tmax_v = col_v_tmax.text_input("Tmax [K]", key="tmax_v")
        ptmax_v = col_v_ptmax.text_input("P(Tmax) [Pa]", key="ptmax_v")

        st.markdown("**Watson (Comp 1):**")
        col_v_hvap, col_v_teb, col_v_tc = st.columns(3)
        hvap_v = col_v_hvap.text_input("Hvap_0 [J/mol]", key="hvap_v")
        teb_v = col_v_teb.text_input("Teb_0 [K]", key="teb_v")
        tc_v = col_v_tc.text_input("Tc [K]", key="tc_v")

        st.markdown("---")     
        st.markdown("**Antoine Estendida (DIPPR 101) (Comp 2):**")    
        col_p_c1, col_p_c2, col_p_c3 = st.columns(3)
        c1_p = col_p_c1.text_input("C1", key="c1_p")
        c2_p = col_p_c2.text_input("C2", key="c2_p")
        c3_p = col_p_c3.text_input("C3", key="c3_p")

        col_p_c4, col_p_c5, col_p_tmin = st.columns(3)
        c4_p = col_p_c4.text_input("C4", key="c4_p")
        c5_p = col_p_c5.text_input("C5", key="c5_p")
        tmin_p = col_p_tmin.text_input("Tmin [K]", key="tmin_p")

        col_p_ptmin, col_p_tmax, col_p_ptmax = st.columns(3)
        ptmin_p = col_p_ptmin.text_input("P(Tmin) [Pa]", key="ptmin_p")
        tmax_p = col_p_tmax.text_input("Tmax [K]", key="tmax_p")
        ptmax_p = col_p_ptmax.text_input("P(Tmax) [Pa]", key="ptmax_p")

        st.markdown('**Watson (Comp 2):**')
        col_p_hvap, col_p_teb, col_p_tc = st.columns(3)
        hvap_p = col_p_hvap.text_input("Hvap_0 [J/mol]", key="hvap_p")
        teb_p = col_p_teb.text_input("Teb_0 [K]", key="teb_p")
        tc_p = col_p_tc.text_input("Tc [K]", key="tc_p")

    with tab2:
        st.markdown("**Coeficientes NRTL:**")
        col_n1, col_n2, col_n3 = st.columns(3)
        a12_nrtl = col_n1.text_input("A12", key="a12_nrtl")
        a21_nrtl = col_n2.text_input("A21", key="a21_nrtl")
        b12_nrtl = col_n3.text_input("B12 [1/K]", key="b12_nrtl")

        col_n4, col_n5, col_n6 = st.columns(3)
        b21_nrtl = col_n4.text_input("B21 [1/K]", key="b21_nrtl")
        c12_nrtl = col_n5.text_input("C12", key="c12_nrtl")
        c21_nrtl = col_n6.text_input("C21", key="c21_nrtl")

        st.markdown("---")
        st.markdown("**Coeficientes Wilson:**")
        col_w1, col_w2, col_w3, col_w4 = st.columns(4)
        a12_wil = col_w1.text_input("A12", key="a12_wil")
        a21_wil = col_w2.text_input("A21", key="a21_wil")
        b12_wil = col_w3.text_input("B12 [1/K]", key="b12_wil")
        b21_wil = col_w4.text_input("B21 [1/K]", key="b21_wil")

    st.markdown("---")
    if st.button("💾 Atribuir Mistura", use_container_width=True):
        v_norm = normalizar_texto(novo_v)
        p_norm = normalizar_texto(novo_p)

        v_clean = mapa_compostos_norm[v_norm][0] if v_norm in mapa_compostos_norm else (novo_v.strip().capitalize() if novo_v else "")
        p_clean = mapa_compostos_norm[p_norm][0] if p_norm in mapa_compostos_norm else (novo_p.strip().capitalize() if novo_p else "")

        if not v_clean or not p_clean:
            st.error("Por favor, preencha os nomes dos dois compostos.")
            return

        c1_v_num, c2_v_num = parse_num(c1_v, 0.0), parse_num(c2_v, 0.0)
        c1_p_num, c2_p_num = parse_num(c1_p, 0.0), parse_num(c2_p, 0.0)
        hvap_v_num, teb_v_num, tc_v_num = parse_num(hvap_v, 0.0), parse_num(teb_v, 0.0), parse_num(tc_v, 0.0)
        hvap_p_num, teb_p_num, tc_p_num = parse_num(hvap_p, 0.0), parse_num(teb_p, 0.0), parse_num(tc_p, 0.0)

        if c1_v_num == 0.0 or c2_v_num == 0.0 or c1_p_num == 0.0 or c2_p_num == 0.0:
            st.error("⚠️ Os parâmetros C1 e C2 de Antoine para ambos os compostos são de preenchimento obrigatório.")
            return

        if hvap_v_num <= 0.0 or teb_v_num <= 0.0 or tc_v_num <= 0.0 or hvap_p_num <= 0.0 or teb_p_num <= 0.0 or tc_p_num <= 0.0:
            st.error("⚠️ Os parâmetros de Watson (Hvap_0, Teb_0 e Tc) para ambos os compostos são de preenchimento obrigatório.")
            return

        if tc_v_num <= teb_v_num or tc_p_num <= teb_p_num:
            st.error("⚠️ A Temperatura Crítica (Tc) deve ser estritamente maior que a Temperatura de Ebulição (Teb_0) para ambos os compostos.")
            return

        nome_mistura_formatado = f"{v_clean}-{p_clean}"

        dict_v_puro = {
            'Antoine': {
                'C1': c1_v_num, 'C2': c2_v_num,
                'C3': parse_num(c3_v, 0.0), 'C4': parse_num(c4_v, 0.0),
                'C5': parse_num(c5_v, 2.0)
            },
            'Limites': {
                'Tmin': parse_num(tmin_v, 150.0), 'P(Tmin)': parse_num(ptmin_v, 1e-3),
                'Tmax': parse_num(tmax_v, 600.0), 'P(Tmax)': parse_num(ptmax_v, 5e6)
            },
            'Watson': {
                'Hvap_0': hvap_v_num,
                'Teb_0': teb_v_num,
                'Tc': tc_v_num
            }
        }

        dict_p_puro = {
            'Antoine': {
                'C1': c1_p_num, 'C2': c2_p_num,
                'C3': parse_num(c3_p, 0.0), 'C4': parse_num(c4_p, 0.0),
                'C5': parse_num(c5_p, 2.0)
            },
            'Limites': {
                'Tmin': parse_num(tmin_p, 150.0), 'P(Tmin)': parse_num(ptmin_p, 1e-3),
                'Tmax': parse_num(tmax_p, 650.0), 'P(Tmax)': parse_num(ptmax_p, 5e6)
            },
            'Watson': {
                'Hvap_0': hvap_p_num,
                'Teb_0': teb_p_num,
                'Tc': tc_p_num
            }
        }

        dict_mistura_bin = {
            'Volatil': v_clean,
            'Pesado': p_clean,
            'NRTL': {
                'a12': parse_num(a12_nrtl, 0.0), 'a21': parse_num(a21_nrtl, 0.0),
                'b12': parse_num(b12_nrtl, 0.0), 'b21': parse_num(b21_nrtl, 0.0),
                'c12': parse_num(c12_nrtl, 0.30), 'c21': parse_num(c21_nrtl, 0.30)
            },
            'Wilson': {
                'a12': parse_num(a12_wil, 0.0), 'a21': parse_num(a21_wil, 0.0),
                'b12': parse_num(b12_wil, 0.0), 'b21': parse_num(b21_wil, 0.0)
            }
        }

        st.session_state['compostos_puros'][v_clean] = dict_v_puro
        st.session_state['compostos_puros'][p_clean] = dict_p_puro
        st.session_state['misturas'][nome_mistura_formatado] = dict_mistura_bin
        st.session_state['mistura_selecionada'] = nome_mistura_formatado

        st.rerun()

#---------------------------------------------- CONFIGURACAO DA BARRA LATERAL -----------------------------------------#
with st.sidebar:
    if st.button("➕ Adicionar Nova Mistura", use_container_width=True):
        modal_adicionar_mistura()

    with st.expander(label='Parâmetros de Operação & Termodinâmica:', expanded=True):
        with st.container(border=True):
            opcoes_misturas = list(MISTURAS.keys())
            if 'mistura_selecionada' not in st.session_state or st.session_state['mistura_selecionada'] not in opcoes_misturas:
                st.session_state['mistura_selecionada'] = opcoes_misturas[0]

            idx_atual = opcoes_misturas.index(st.session_state['mistura_selecionada'])

            Mixture = st.selectbox(
                label="Mistura:",
                options=opcoes_misturas,
                index=idx_atual
            )
            st.session_state['mistura_selecionada'] = Mixture

            Nome_Volatil = MISTURAS[Mixture]['Volatil']
            Nome_Pesado = MISTURAS[Mixture]['Pesado']
            Comp_Volatil = COMPOSTOS_PUROS[Nome_Volatil]
            Comp_Pesado = COMPOSTOS_PUROS[Nome_Pesado]

            Thermo_model = st.selectbox(label='Modelo termodinâmico:', options=db.MODELOS)

            P_min_Pa = max(Comp_Volatil['Limites']['P(Tmin)'], Comp_Pesado['Limites']['P(Tmin)'])
            P_min_kPa = float(np.ceil(P_min_Pa / 1000.0 * 100) / 100)

            Pressao_kPa = st.slider(
                label='Pressão de Operação da Coluna: [kPa]',
                min_value=P_min_kPa,
                max_value=506.625,
                value=max(101.325, P_min_kPa),
                step=1.0,
                help='Pressão total mantida constante na coluna de destilação.'
            )
            P_total_Pa = Pressao_kPa * 1000.0
            st.caption(f"Equivalente: **{Pressao_kPa / 101.325:.2f} atm**")

        with st.container(border=True):
            st.markdown("**Análise Isotérmica (P-x-y):**")
            T_inf = max(Comp_Volatil['Limites']['Tmin'], Comp_Pesado['Limites']['Tmin'])
            T_sup = min(Comp_Pesado['Limites']['Tmax'], Comp_Volatil['Limites']['Tmax'])

            Temperatura = st.slider(
                label='Temperatura de Referência (P-x-y): [K]',
                min_value=T_inf,
                max_value=T_sup,
                value=(T_inf + T_sup) / 2,
                step=0.02
            )

    x_azeo, tipo_azeo = be.detectar_azeotropo_isobarico(
        P_total_Pa, Thermo_model, MISTURAS[Mixture], COMPOSTOS_PUROS, Nome_Volatil, Nome_Pesado
    )

    with st.expander(label='McCabe-Thiele:', expanded=True):
        max_xd_permitido = 0.99
        min_xb_permitido = 0.01

        if x_azeo is not None:
            if tipo_azeo == 'minimo':
                max_xd_permitido = max(0.51, min(0.99, round(x_azeo - 0.01, 2)))
                st.warning(
                    f"⚠️ Azeótropo de T_eb mínimo em x = {x_azeo:.3f}! "
                    f"Topo limitado a Xd ≤ {max_xd_permitido:.2f}."
                )
            else:
                min_xb_permitido = min(0.95, max(0.01, round(x_azeo + 0.01, 2)))
                st.warning(
                    f"⚠️ Azeótropo de T_eb máximo em x = {x_azeo:.3f}! "
                    f"Fundo limitado a Xb ≥ {min_xb_permitido:.2f}."
                )

        with st.container(border=True):
            Vazao_F = st.number_input(
                label='Vazão de Alimentação (F) [kmol/h]:',
                min_value=1.0, 
                max_value=10000.0, 
                value=100.0, 
                step=1.0
            )
            min_xd_slider = max(0.50, round(min_xb_permitido + 0.03, 2))
            Comp_Topo = st.slider(
                label=f'Fração molar de {Nome_Volatil} no topo (Xd)',
                min_value=min_xd_slider, 
                max_value=max_xd_permitido, 
                value=min(max(0.90, min_xd_slider), max_xd_permitido), 
                step=0.01
            )
            min_zf_slider = max(0.05, round(min_xb_permitido + 0.01, 2))
            max_zf_slider = round(Comp_Topo - 0.01, 2)
            Comp_Alimentacao = st.slider(
                label=f'Fração molar de {Nome_Volatil} na alimentação (Zf)',
                min_value=min_zf_slider, 
                max_value=max_zf_slider, 
                value=min(max(0.50, min_zf_slider + 0.01), max_zf_slider), 
                step=0.01
            )
            max_xb_slider = round(Comp_Alimentacao - 0.01, 2)
            Comp_Fundo = st.slider(
                label=f'Fração molar de {Nome_Volatil} no fundo (Xb)',
                min_value=min_xb_permitido, 
                max_value=max_xb_slider, 
                value=min(max(0.05, min_xb_permitido), max_xb_slider), 
                step=0.01
            )

        with st.container(border=True):
            Refluxo = st.slider(
                label='Razão de Refluxo (R)', 
                min_value=0.00, 
                max_value=10.00, 
                value=1.00, 
                step=0.01
            )
            Parametro_Q = st.slider(
                label='Parâmetro q', 
                min_value=0.00, 
                max_value=1.00, 
                value=1.00, 
                step=0.01,
                help='1.0 = Líquido Saturado | 0.0 = Vapor Saturado'
            )
            Tipo_Murphree = st.radio(
                label="Base da Eficiência de Murphree:",
                options=['Vapor (EMV)', 'Líquido (EML)'],
                index=0,
                horizontal=True,
                help="'Vapor' é o padrão industrial mais utilizado. 'Líquido' aplica a eficiência no rebate horizontal."
            )
            Murphree = st.slider(
                label=f'Eficiência ({Tipo_Murphree})', min_value=0.01, max_value=1.00, value=1.00, step=0.01
            )

    
    with st.container(border=True):              
        st.markdown("**Análise de Sensibilidade & Relatório**")
        
        rodar_sensibilidade = st.button(
            "Gerar Análise de Sensibilidade",
            use_container_width=True,
            help="Clique para processar as varreduras estáticas de Refluxo, Pressão e Parâmetro q."
        )
                       
           

#---------------------------------------------- PROCESSAMENTO DOS CALCULOS PRINCIPAIS ---------------------------------#
x1_grid, y1_isobarico, T_bolha_arr = be.gerar_diagrama_TXY(
    P_total_Pa=P_total_Pa,
    modelo=Thermo_model,
    ctes_thermo=MISTURAS[Mixture],
    db_compostos=COMPOSTOS_PUROS,
    nome_v=Nome_Volatil,
    nome_p=Nome_Pesado
)

mc_dados = be.calcular_mccabe_thiele_isobarico(
    xd=Comp_Topo, 
    xb=Comp_Fundo, 
    zf=Comp_Alimentacao,
    R=Refluxo, 
    q=Parametro_Q, 
    murphree=Murphree, 
    tipo_murphree=Tipo_Murphree,
    P_total_Pa=P_total_Pa, 
    modelo=Thermo_model, 
    ctes_thermo=MISTURAS[Mixture],
    db_compostos=COMPOSTOS_PUROS, 
    nome_v=Nome_Volatil, 
    nome_p=Nome_Pesado
)

balanco = be.calcular_balanco_termico(
    F_kmol_h=Vazao_F, 
    xd=Comp_Topo, 
    xb=Comp_Fundo, 
    zf=Comp_Alimentacao,
    R=Refluxo, 
    q=Parametro_Q, 
    P_total_Pa=P_total_Pa, 
    modelo=Thermo_model,
    ctes_thermo=MISTURAS[Mixture], 
    db_compostos=COMPOSTOS_PUROS,
    nome_v=Nome_Volatil, 
    nome_p=Nome_Pesado
)

R_min = be.calcular_refluxo_minimo_isobarico(
    xd=Comp_Topo,
    zf=Comp_Alimentacao,
    q=Parametro_Q,
    P_total_Pa=P_total_Pa,
    modelo=Thermo_model,
    ctes_thermo=MISTURAS[Mixture],
    db_compostos=COMPOSTOS_PUROS,
    nome_v=Nome_Volatil,
    nome_p=Nome_Pesado
)

n_min = be.calcular_estagios_minimos_isobarico(
    xd=Comp_Topo,
    xb=Comp_Fundo,
    P_total_Pa=P_total_Pa,
    modelo=Thermo_model,
    ctes_thermo=MISTURAS[Mixture],
    db_compostos=COMPOSTOS_PUROS,
    nome_v=Nome_Volatil,
    nome_p=Nome_Pesado
)

#---------------------------------------------- RENDERIZACAO MCCABE-THIELE --------------------------------------------#
Curva_equilibrio = go.Figure()

Curva_equilibrio.add_trace(go.Scatter(
    x=x1_grid, 
    y=y1_isobarico, 
    mode='lines',
    name=f'Equilíbrio ({Pressao_kPa:.1f} kPa)',
    hovertemplate='<b>X:</b> %{x:.3f}<br><b>Y:</b> %{y:.3f}<extra></extra>',
    line=dict(color='#E60505', width=2)
))


Curva_equilibrio.add_trace(go.Scatter(
    x=[0, 1], 
    y=[0, 1], 
    mode='lines', 
    name='Linha 45°',
    line=dict(color='gray', dash='dash')
))
Curva_equilibrio.add_trace(go.Scatter(
    x=mc_dados['x_q'], 
    y=mc_dados['y_q'], 
    mode='lines', 
    name='Linha q',
    hovertemplate='<b>X:</b> %{x:.3f}<br><b>Y:</b> %{y:.3f}<extra></extra>',
    line=dict(color='green', width=2)
))
Curva_equilibrio.add_trace(go.Scatter(
    x=mc_dados['x_losr'], 
    y=mc_dados['y_losr'], 
    mode='lines', 
    name='Linha de retificação',
    line=dict(color='orange', width=2)
))
Curva_equilibrio.add_trace(go.Scatter(
    x=mc_dados['x_lose'], 
    y=mc_dados['y_lose'], 
    mode='lines', 
    name='Linha esgotamento',
    line=dict(color='purple', width=2)
))
Curva_equilibrio.add_trace(go.Scatter(
    x=mc_dados['x_degraus'], 
    y=mc_dados['y_degraus'], 
    mode='lines', 
    name='Estágios',
    hovertemplate='<b>x:</b> %{x:.3f}<br><b>y:</b> %{y:.3f}<extra></extra>',
    line=dict(color='#0055FF', width=1.5)
))

Curva_equilibrio.update_layout(
    title=f"Diagrama McCabe-Thiele Isobárico ({Pressao_kPa:.1f} kPa) | Estágios Calculados: {mc_dados['n_estagios']}",
    xaxis_title=f'<b>Fração líquida {Nome_Volatil}</b>',
    yaxis_title=f'<b>Fração gasosa {Nome_Volatil}</b>',
    height=670, margin=dict(l=20, r=20, t=40, b=20),
    hovermode='x unified',
    xaxis=dict(range=[0, 1], dtick=0.05, showgrid=True, gridcolor='#5C768D'),
    yaxis=dict(range=[0, 1], dtick=0.05, showgrid=True, gridcolor='#5C768D')
)

st.plotly_chart(Curva_equilibrio, use_container_width=True)

if mc_dados['nao_convergiu']:
    st.error("⚠️ Os estágios não atingiram a composição de fundo Xb. Aumente a Razão de Refluxo para operar acima do ponto de Pinch.")

#---------------------------------------------- PAINEL DE METRICAS OPERACIONAIS ----------------------------------------#

st.markdown('-----')
col_1, col_2, col_3, col_4, col_5 = st.columns(5)
with col_1:
    est_alim = mc_dados.get('estagio_alimentacao')
    st.metric(
        "N° de Estágios Reais",
        f"{mc_dados['n_estagios']}",
        delta=f"Alim.: {est_alim}º estágio (topo → fundo)" if est_alim is not None else "Alim.: N/A",
        delta_color="off"
    )
    st.metric("Vazão Destilado (D)", f"{balanco['D']:.2f} kmol/h")
    st.caption(f"Temp. Topo: **{balanco['T_topo_C']:.1f} °C**")
with col_2:
    st.metric("Eficiência Apl.", f"{Murphree*100:.0f}%")
    st.metric("Vazão Resíduo (B)", f"{balanco['B']:.2f} kmol/h")
    st.caption(f"Temp. Fundo: **{balanco['T_fundo_C']:.1f} °C**")
with col_3:
    st.metric("Pressão da Coluna", f"{Pressao_kPa:.2f} kPa")
    st.metric("Condensador (Qc)", f"{balanco['Qc_kW']:.2f} kW")
    st.caption(f"ΔHvap topo: **{balanco['hvap_topo_kJ_mol']:.2f} kJ/mol**")
with col_4:
    st.metric(f"Temp. Ebulição {Nome_Volatil} (Puro)", f"{T_bolha_arr[-1] - 273.15:.1f} °C")
    st.metric("Refervedor (Qr)", f"{balanco['Qr_kW']:.2f} kW")
    st.caption(f"ΔHvap fundo: **{balanco['hvap_fundo_kJ_mol']:.2f} kJ/mol**")
with col_5:
    r_min_display = f"{R_min:.4f}" if np.isfinite(R_min) else "Pinch"
    st.metric("Razão Refluxo Mínima.", r_min_display)
    st.metric("Número mínimo de estágios.", f"{n_min}")
    if np.isfinite(R_min) and R_min > 0:
        st.caption(f"Relação R / Rmín: **{Refluxo / R_min:.2f}x**")

#---------------------------------------------- RENDERIZACAO DIAGRAMA T-X-Y -------------------------------------------#
Diagrama_TXY = go.Figure()
T_bolha_C = T_bolha_arr - 273.15

Diagrama_TXY.add_trace(go.Scatter(
    x=x1_grid, 
    y=T_bolha_C, 
    mode='lines', 
    name='Ponto de Bolha',
    hovertemplate='<b>x:</b> %{x:.3f}<br><b>T Bolha:</b> %{y:.2f} °C<extra></extra>',
    line=dict(color='blue', width=2)
))
Diagrama_TXY.add_trace(go.Scatter(
    x=y1_isobarico, 
    y=T_bolha_C, 
    mode='lines', 
    name='Ponto de Orvalho',
    hovertemplate='<b>y:</b> %{x:.3f}<br><b>T Orvalho:</b> %{y:.2f} °C<extra></extra>',
    line=dict(color='#E60505', width=2)
))
Diagrama_TXY.update_layout(
    title=f'Diagrama T-x-y Isobárico ({Pressao_kPa:.2f} kPa)',
    xaxis_title=f'<b>Composição X, Y ({Nome_Volatil})</b>',
    yaxis_title='<b>Temperatura [°C]</b>',
    hovermode='x unified', height=670,
    margin=dict(l=20, r=20, t=40, b=20),
    xaxis=dict(range=[0, 1], dtick=0.05, showgrid=True, gridcolor='#5C768D'),
    yaxis=dict(showgrid=True, gridcolor='#5C768D')
)
st.plotly_chart(Diagrama_TXY, use_container_width=True)

#---------------------------------------------- RENDERIZACAO DIAGRAMA P-X-Y -------------------------------------------#
Psat_V_iso = to.Psat(T_K=Temperatura, **Comp_Volatil['Antoine'])
Psat_P_iso = to.Psat(T_K=Temperatura, **Comp_Pesado['Antoine'])
Gamma_V_iso, Gamma_P_iso = to.Coef_atv(Thermo_model, MISTURAS[Mixture], Temperatura)
Y_v_iso, Ptotal_iso = be.Raoult_mod(x1_grid, Gamma_V_iso, Gamma_P_iso, Psat_V_iso, Psat_P_iso)

Diagrama_PXY = go.Figure()
Diagrama_PXY.add_trace(go.Scatter(
    x=x1_grid, 
    y=Ptotal_iso / 1000, 
    mode='lines', 
    name='Ponto de Bolha',
    line=dict(color='blue')
))
Diagrama_PXY.add_trace(go.Scatter(
    x=Y_v_iso, 
    y=Ptotal_iso / 1000, 
    mode='lines', 
    name='Ponto de Orvalho',
    line=dict(color='#E60505')
))
Diagrama_PXY.update_layout(
    title=f'Diagrama P-x-y Isotérmico (Ref: {Temperatura:.2f} K / {Temperatura - 273.15:.2f} °C)',
    xaxis_title=f'<b>Composição X, Y ({Nome_Volatil})</b>',
    yaxis_title='<b>Pressão Total [kPa]</b>',
    hovermode='x unified', height=670,
    margin=dict(l=20, r=20, t=40, b=20),
    xaxis=dict(range=[0, 1], dtick=0.05, showgrid=True, gridcolor='#5C768D'),
    yaxis=dict(showgrid=True, gridcolor='#5C768D')
)
st.plotly_chart(Diagrama_PXY, use_container_width=True)

#---------------------------------------------- RENDERIZACAO ANALISE SENSIBILIDADE ------------------------------------#
if rodar_sensibilidade:
    with st.spinner("⏳ Processando análise de sensibilidade..."):
        st.session_state['dados_sensibilidade'] = be.gerar_analise_sensibilidade(
            xd=Comp_Topo, 
            xb=Comp_Fundo, 
            zf=Comp_Alimentacao,
            R_ref=Refluxo, 
            q_ref=Parametro_Q, 
            murphree=Murphree,
            tipo_murphree=Tipo_Murphree, 
            P_total_Pa_atual=P_total_Pa,
            modelo=Thermo_model, 
            ctes_thermo=MISTURAS[Mixture],
            db_compostos=COMPOSTOS_PUROS, 
            nome_v=Nome_Volatil, 
            nome_p=Nome_Pesado,
            F_kmol_h=Vazao_F
        )

if 'dados_sensibilidade' in st.session_state and st.session_state['dados_sensibilidade'] is not None:
    sensib = st.session_state['dados_sensibilidade']

    fig1 = make_subplots(specs=[[{"secondary_y": True}]])
    fig1.add_trace(
        go.Scatter(
            x=sensib['R_vec'], 
            y=sensib['n_estagios_R'],
            name="Nº de Estágios Reais", 
            mode='lines+markers',
            line=dict(color='royalblue', width=2.5, shape='hv')
        ),
        secondary_y=False
    )
    fig1.add_trace(
        go.Scatter(
            x=sensib['R_vec'], 
            y=sensib['Qc_kW_R'],
            name="Carga Térmica no Condensador ", 
            mode='lines+markers',
            line=dict(color='purple', width=2.5, dash='dash')
        ),
        secondary_y=True
    )
    fig1.add_vline(
        x=Refluxo, 
        line_width=1.5, 
        line_dash="dot", 
        line_color="gray",
        annotation_text=f"R Atual ({Refluxo:.2f})"
    )
    fig1.update_layout(
        title="1. Otimização: Razão de Refluxo vs. Estágios & Carga Térmica no Condensador",
        xaxis_title="Razão de Refluxo", 
        template="plotly_white", 
        height=600,
        hovermode="x unified", 
        legend=dict(x=0.01, y=1.12, orientation="h")
    )
    fig1.update_yaxes(
        title_text="<b>Nº de Estágios Reais</b>", 
        title_font=dict(color="royalblue"), 
        secondary_y=False
        )
    fig1.update_yaxes(
        title_text="<b>Carga Térmica Qc (kW)</b>", 
        title_font=dict(color="purple"), 
        secondary_y=True
        )

    st.plotly_chart(fig1, use_container_width=True)

    fig2 = make_subplots(specs=[[{"secondary_y": True}]])
    fig2.add_trace(
        go.Scatter(
            x=sensib['P_vec_kPa'], 
            y=sensib['T_topo_C_P'],
            name="Temp. Topo", 
            mode='lines+markers',
            line=dict(color='deepskyblue', width=2)
        ),
        secondary_y=False
    )
    fig2.add_trace(
        go.Scatter(
            x=sensib['P_vec_kPa'], 
            y=sensib['T_fundo_C_P'],
            name="Temp. Fundo", 
            mode='lines+markers',
            line=dict(color='orange', width=2)
        ),
        secondary_y=False
    )
    fig2.add_trace(
        go.Scatter(
            x=sensib['P_vec_kPa'], 
            y=sensib['Qc_kW_P'],
            name="Carga Térmica no Condensador", 
            mode='lines+markers',
            line=dict(color='purple', width=2.5, dash='dash')
        ),
        secondary_y=True
    )
    fig2.add_vline(
        x=Pressao_kPa, 
        line_width=1.5, 
        line_dash="dot", 
        line_color="gray",
        annotation_text=f"P Atual ({Pressao_kPa:.1f} kPa)"
    )
    fig2.update_layout(
        title="2. Efeito da Pressão de Operação nas Temperaturas e Carga Térmica no Condensador",
        xaxis_title="Pressão da Coluna (kPa)", 
        template="plotly_white", 
        height=600,
        hovermode="x unified", 
        legend=dict(x=0.01, y=1.12, orientation="h")
    )
    fig2.update_yaxes(
        title_text="<b>Temperatura (°C)</b>", 
        secondary_y=False
        )
    fig2.update_yaxes(
        title_text="<b>Carga Térmica Qc (kW)</b>", 
        title_font=dict(color="purple"),
        secondary_y=True
        )

    st.plotly_chart(fig2, use_container_width=True)

    fig3 = make_subplots(specs=[[{"secondary_y": True}]])
    fig3.add_trace(
        go.Scatter(
            x=sensib['q_vec'], 
            y=sensib['Qr_kW_q'],
            name="Carga Térmica Refervedor", 
            mode='lines+markers',
            line=dict(color='crimson', width=2.5)
        ),
        secondary_y=False
    )
    fig3.add_trace(
        go.Scatter(
            x=sensib['q_vec'], 
            y=sensib['n_estagios_q'],
            name="Nº de Estágios Reais", 
            mode='lines+markers',
            line=dict(color='royalblue', width=2, dash='dot', shape='hv')
        ),
        secondary_y=True
    )
    fig3.add_vline(
        x=Parametro_Q, 
        line_width=1.5, 
        line_dash="dot", 
        line_color="gray",
        annotation_text=f"q Atual ({Parametro_Q:.2f})"
    )
    fig3.update_layout(
        title="3. Condição Térmica da Alimentação vs. Carga térmica no Refervedor & Estágios",
        xaxis_title="Parâmetro q (0 = Vapor Sat. | 1 = Líquido Sat.)", 
        template="plotly_white", 
        height=600,
        hovermode="x unified", legend=dict(x=0.01, y=1.12, orientation="h")
    )
    fig3.update_yaxes(
        title_text="<b>Carga Térmica Qr (kW)</b>", 
        title_font=dict(color="crimson"), 
        secondary_y=False
        )
    fig3.update_yaxes(
        title_text="<b>Nº de Estágios Reais</b>", 
        title_font=dict(color="royalblue"), 
        secondary_y=True, 
        tickformat='d'
        )

    st.plotly_chart(fig3, use_container_width=True)

#---------------------------------------------- EXPORTACAO DO RELATORIO E DADOS (.ZIP) ---------------------------------#
sensib_data_pdf = st.session_state.get('dados_sensibilidade', None)
fig1_pdf = fig1 if 'fig1' in locals() else None
fig2_pdf = fig2 if 'fig2' in locals() else None
fig3_pdf = fig3 if 'fig3' in locals() else None

pdf_bytes = be.gerar_pdf_relatorio(
    mistura=Mixture, 
    modelo=Thermo_model, 
    P_kPa=Pressao_kPa, 
    F=Vazao_F, 
    xd=Comp_Topo, 
    zf=Comp_Alimentacao, 
    xb=Comp_Fundo, 
    R=Refluxo, 
    q=Parametro_Q, 
    murphree=Murphree, 
    tipo_murphree=Tipo_Murphree, 
    mc_dados=mc_dados, 
    balanco=balanco, 
    fig_mc=Curva_equilibrio, 
    fig_txy=Diagrama_TXY, 
    fig_pxy=Diagrama_PXY, 
    x1_grid=x1_grid, 
    y1_isobarico=y1_isobarico, 
    T_bolha_arr=T_bolha_arr, 
    Y_v_iso=Y_v_iso, 
    Ptotal_iso=Ptotal_iso,
    T_ref_K=Temperatura,  
    sensib_dados=sensib_data_pdf,
    fig_sensib_R=fig1_pdf,
    fig_sensib_P=fig2_pdf,
    fig_sensib_q=fig3_pdf,
    R_min=R_min,
    n_min=n_min
)

zip_bytes = be.gerar_pacote_zip(
    pdf_bytes=pdf_bytes,
    mistura=Mixture,
    modelo=Thermo_model,
    P_kPa=Pressao_kPa,
    F=Vazao_F,
    xd=Comp_Topo,
    zf=Comp_Alimentacao,
    xb=Comp_Fundo,
    R=Refluxo,
    q=Parametro_Q,
    murphree=Murphree,
    tipo_murphree=Tipo_Murphree,
    mc_dados=mc_dados,
    balanco=balanco,
    x1_grid=x1_grid,
    y1_isobarico=y1_isobarico,
    T_bolha_arr=T_bolha_arr,
    Y_v_iso=Y_v_iso,
    Ptotal_iso=Ptotal_iso,
    T_ref_K=Temperatura,
    sensib_dados=sensib_data_pdf,
    R_min=R_min,
    n_min=n_min
)

caminho_manual = os.path.join(diretorio_atual, "Manual_YuNno_thermo.pdf")

@st.cache_data(show_spinner=False)
def carregar_manual_pdf(caminho):
    if os.path.exists(caminho):
        with open(caminho, "rb") as f:
            return f.read()
    return None

manual_bytes = carregar_manual_pdf(caminho_manual)

with st.sidebar:
    st.download_button(
        label="Baixar relatório e dados (.ZIP)",
        data=zip_bytes,
        file_name=f"YuNno_Thermo_{Mixture}_{Pressao_kPa:.0f}kPa.zip",
        mime="application/zip",
        use_container_width=True
    )
    if manual_bytes is not None:
        st.download_button(
            label="Manual do Usuário (.PDF)",
            data=manual_bytes,
            file_name="Manual_YuNno_thermo.pdf",
            mime="application/pdf",
            use_container_width=True
        )
    else:
        st.button(
            label="Manual do Usuário (.PDF)",
            disabled=True,
            use_container_width=True,
            help="Arquivo Manual_Usuario_YuNno_Thermo.pdf não encontrado no diretório da aplicação."
        )