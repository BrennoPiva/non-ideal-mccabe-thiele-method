#---------------------------------------------------- MOTOR DE CALCULOS E PDF ------------------------------------------#
import tools as to
import numpy as np
import scipy.optimize as sp_opt
import streamlit as st
import io
import copy
from reportlab.lib.pagesizes import letter, landscape
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image as RLImage, KeepTogether, PageBreak
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib import colors
import io
import copy
import zipfile
import os


#------------------------------------------------- EQUILIBRIO LIQUIDO-VAPOR P-X-Y --------------------------------------#
def Raoult_mod(x1, gamma1, gamma2, Psat_A, Psat_B):
   
    x2 = 1.0 - x1
    P1_parcial = x1 * gamma1 * Psat_A
    P2_parcial = x2 * gamma2 * Psat_B
    P = P1_parcial + P2_parcial
    y1 = P1_parcial / P
    return y1, P

#---------------------------------------------------- EQUILIBRIO ISOBARICO T-X-Y ---------------------------------------#
def obter_y_equilibrio_isobarico(x_val, P_total_Pa, modelo, ctes_thermo, db_compostos, nome_v, nome_p):
   
    _, y1 = to.temperatura_bolha(x_val, P_total_Pa, modelo, ctes_thermo, db_compostos, nome_v, nome_p)
    return y1

def obter_x_equilibrio_isobarico(y_target, P_total_Pa, modelo, ctes_thermo, db_compostos, nome_v, nome_p):
   
    def erro_eq(x):
        return obter_y_equilibrio_isobarico(x, P_total_Pa, modelo, ctes_thermo, db_compostos, nome_v, nome_p) - y_target

    try:
        x_sol = sp_opt.brentq(erro_eq, 1e-6, 1.0 - 1e-6)
        return float(x_sol)
    except ValueError:
        return float(y_target)

def detectar_azeotropo_isobarico(P_total_Pa, modelo, ctes_thermo, db_compostos, nome_v, nome_p):
   
    def erro_azeo(x):
        y = obter_y_equilibrio_isobarico(x, P_total_Pa, modelo, ctes_thermo, db_compostos, nome_v, nome_p)
        return y - x

    try:
        x_azeo = float(sp_opt.brentq(erro_azeo, 0.02, 0.98))
        
        x_teste = max(0.01, x_azeo - 0.02)
        y_teste = obter_y_equilibrio_isobarico(x_teste, P_total_Pa, modelo, ctes_thermo, db_compostos, nome_v, nome_p)
        tipo_azeo = 'minimo' if y_teste > x_teste else 'maximo'
        return x_azeo, tipo_azeo
    except ValueError:
        return None, None

def gerar_diagrama_TXY(P_total_Pa, modelo, ctes_thermo, db_compostos, nome_v, nome_p):
   
    x_grid = np.linspace(0.00, 1.00, 101)
    T_bolha_list = []
    y_list = []

    for x_val in x_grid:
        T_b, y_b = to.temperatura_bolha(x_val, P_total_Pa, modelo, ctes_thermo, db_compostos, nome_v, nome_p)
        T_bolha_list.append(T_b)
        y_list.append(y_b)

    return x_grid, np.array(y_list), np.array(T_bolha_list)

#------------------------------------------- LIMITES TERMODINAMICOS (Rmin e Nmin) --------------------------------------#
def calcular_refluxo_minimo_isobarico(xd, zf, q, P_total_Pa, modelo, ctes_thermo, db_compostos, nome_v, nome_p):
    
    if abs(q - 1.0) < 1e-4:
        xq_star = zf
        yq_star = obter_y_equilibrio_isobarico(xq_star, P_total_Pa, modelo, ctes_thermo, db_compostos, nome_v, nome_p)
    else:
        a2 = q / (q - 1.0)
        b2 = -zf / (q - 1.0)
        def erro_intersec_q(x):
            y_eq = obter_y_equilibrio_isobarico(x, P_total_Pa, modelo, ctes_thermo, db_compostos, nome_v, nome_p)
            return y_eq - (a2 * x + b2)

        x_busca_min = 1e-4
        x_busca_max = max(1e-3, xd - 1e-4)
        try:
            xq_star = float(sp_opt.brentq(erro_intersec_q, x_busca_min, x_busca_max))
        except ValueError:
            xq_star = zf
        yq_star = obter_y_equilibrio_isobarico(xq_star, P_total_Pa, modelo, ctes_thermo, db_compostos, nome_v, nome_p)


    inclinacao_max = (xd - yq_star) / (xd - xq_star) if abs(xd - xq_star) > 1e-6 else 0.0

   
    if xd - xq_star > 0.02:
        x_tang_grid = np.linspace(xq_star, xd - 0.005, 40)
        for xt in x_tang_grid:
            yt = obter_y_equilibrio_isobarico(xt, P_total_Pa, modelo, ctes_thermo, db_compostos, nome_v, nome_p)
            inclinacao_cand = (xd - yt) / (xd - xt)
            if inclinacao_cand > inclinacao_max:
                inclinacao_max = inclinacao_cand

    if inclinacao_max >= 1.0 - 1e-6:
        return float('inf')

    R_min = inclinacao_max / (1.0 - inclinacao_max)
    return float(max(0.0, R_min))


def calcular_estagios_minimos_isobarico(xd, xb, P_total_Pa, modelo, ctes_thermo, db_compostos, nome_v, nome_p):
    
    x_atual = xd
    y_atual = xd
    n_min = 0
    max_iter = 100

    while x_atual > xb and n_min < max_iter:
        n_min += 1
        x_eq = obter_x_equilibrio_isobarico(y_atual, P_total_Pa, modelo, ctes_thermo, db_compostos, nome_v, nome_p)
        if x_eq >= x_atual - 1e-5:
            break
        x_atual = x_eq
        y_atual = x_eq

    return int(n_min)

#---------------------------------------------- DIMENSIONAMENTO MCCABE-THIELE -----------------------------------------#
def calcular_mccabe_thiele_isobarico(xd, xb, zf, R, q, murphree, tipo_murphree, P_total_Pa, modelo, ctes_thermo, db_compostos, nome_v, nome_p):
    if abs(q - 1.0) < 1e-4:
        x_q = zf
        y_q = (R / (R + 1.0)) * x_q + (xd / (R + 1.0))
    else:
        a1 = R / (R + 1.0)
        b1 = xd / (R + 1.0)
        a2 = q / (q - 1.0)
        b2 = -zf / (q - 1.0)
        x_q = (b2 - b1) / (a1 - a2)
        y_q = a1 * x_q + b1

    m_lose = (y_q - xb) / (x_q - xb) if abs(x_q - xb) > 1e-6 else 1.0

    def y_losr(x):
        return (R / (R + 1.0)) * x + (xd / (R + 1.0))

    def y_lose(x):
        return m_lose * (x - xb) + xb

    x_degraus = [xd]
    y_degraus = [xd]

    x_atual = xd
    y_atual = xd
    n_estagios = 0
    estagio_alimentacao = None
    max_iter = 100
    nao_convergiu = False

    while x_atual > xb and n_estagios < max_iter:
        n_estagios += 1
        x_eq_ideal = obter_x_equilibrio_isobarico(y_atual, P_total_Pa, modelo, ctes_thermo, db_compostos, nome_v, nome_p)

        if tipo_murphree == 'Líquido (EML)':
            x_prato = x_atual - murphree * (x_atual - x_eq_ideal)
            x_degraus.append(x_prato)
            y_degraus.append(y_atual)

            if estagio_alimentacao is None and x_prato <= x_q:
                estagio_alimentacao = n_estagios

            if x_prato <= xb:
                x_degraus.append(x_prato)
                y_degraus.append(x_prato)
                break

            y_op = y_losr(x_prato) if (x_prato >= x_q) else y_lose(x_prato)
            x_degraus.append(x_prato)
            y_degraus.append(y_op)

            if x_prato >= x_atual - 1e-5:
                nao_convergiu = True
                break

            x_atual = x_prato
            y_atual = y_op

        else:
            x_degraus.append(x_eq_ideal)
            y_degraus.append(y_atual)

            if estagio_alimentacao is None and x_eq_ideal <= x_q:
                estagio_alimentacao = n_estagios

            if x_eq_ideal <= xb:
                x_degraus.append(x_eq_ideal)
                y_degraus.append(x_eq_ideal)
                break

            y_op_ideal = y_losr(x_eq_ideal) if (x_eq_ideal >= x_q) else y_lose(x_eq_ideal)
            y_real = y_atual - murphree * (y_atual - y_op_ideal)

            x_degraus.append(x_eq_ideal)
            y_degraus.append(y_real)

            if x_eq_ideal >= x_atual - 1e-5:
                nao_convergiu = True
                break

            x_atual = x_eq_ideal
            y_atual = y_real

    if x_atual > xb and n_estagios >= max_iter:
        nao_convergiu = True

    if nao_convergiu:
        estagio_alimentacao = None
    elif estagio_alimentacao is None and n_estagios > 0:
        estagio_alimentacao = n_estagios

    x_losr_plot = np.array([x_q, xd])
    y_losr_plot = y_losr(x_losr_plot)

    x_lose_plot = np.array([xb, x_q])
    y_lose_plot = y_lose(x_lose_plot)

    x_q_plot = np.array([zf, x_q])
    y_q_plot = np.array([zf, y_q])         

    return {
        'x_degraus': x_degraus,
        'y_degraus': y_degraus,
        'x_losr': x_losr_plot,
        'y_losr': y_losr_plot,
        'x_lose': x_lose_plot,
        'y_lose': y_lose_plot,
        'x_q': x_q_plot,
        'y_q': y_q_plot,       
        'n_estagios': n_estagios,
        'estagio_alimentacao': estagio_alimentacao,
        'nao_convergiu': nao_convergiu
    }

#---------------------------------------------- BALANCO DE MASSA E ENERGIA ---------------------------------------------#
def calcular_balanco_termico(F_kmol_h, xd, xb, zf, R, q, P_total_Pa, modelo, ctes_thermo, db_compostos, nome_v, nome_p):
    
    D = F_kmol_h * (zf - xb) / (xd - xb)
    B = F_kmol_h - D

    V_topo = (R + 1.0) * D
    V_fundo = V_topo - (1.0 - q) * F_kmol_h

    T_topo, _ = to.temperatura_bolha(xd, P_total_Pa, modelo, ctes_thermo, db_compostos, nome_v, nome_p)
    T_fundo, _ = to.temperatura_bolha(xb, P_total_Pa, modelo, ctes_thermo, db_compostos, nome_v, nome_p)

    hvap_topo_J_mol = to.hvap_mistura_watson(T_topo, xd, db_compostos, nome_v, nome_p)
    hvap_fundo_J_mol = to.hvap_mistura_watson(T_fundo, xb, db_compostos, nome_v, nome_p)

    Qc_kW = (V_topo * hvap_topo_J_mol) / 3600.0
    Qr_kW = (V_fundo * hvap_fundo_J_mol) / 3600.0

    return {
        'D': D,
        'B': B,
        'V_topo': V_topo,
        'V_fundo': V_fundo,
        'T_topo_C': T_topo - 273.15,
        'T_fundo_C': T_fundo - 273.15,
        'hvap_topo_kJ_mol': hvap_topo_J_mol / 1000.0,
        'hvap_fundo_kJ_mol': hvap_fundo_J_mol / 1000.0,
        'Qc_kW': Qc_kW,
        'Qr_kW': Qr_kW
    }

#---------------------------------------------- VARREDURA DE SENSIBILIDADE --------------------------------------------#
@st.cache_data(show_spinner=False)
def gerar_analise_sensibilidade(xd, xb, zf, R_ref, q_ref, murphree, tipo_murphree, P_total_Pa_atual, modelo, ctes_thermo, db_compostos, nome_v, nome_p, F_kmol_h):
    
    R_vec = np.linspace(0.15, 10.0, 100)
    n_estagios_R = []
    Qc_kW_R = []
    
    for r in R_vec:
        res_mc = calcular_mccabe_thiele_isobarico(xd, xb, zf, r, q_ref, murphree, tipo_murphree, P_total_Pa_atual, modelo, ctes_thermo, db_compostos, nome_v, nome_p)
        res_term = calcular_balanco_termico(F_kmol_h, xd, xb, zf, r, q_ref, P_total_Pa_atual, modelo, ctes_thermo, db_compostos, nome_v, nome_p)
        n_estagios_R.append(res_mc['n_estagios'] if not res_mc['nao_convergiu'] else None)
        Qc_kW_R.append(res_term['Qc_kW'])

    P_vec_kPa = np.linspace(10.0, 600.0, 100)
    T_topo_C_P, T_fundo_C_P, Qc_kW_P = [], [], []
    
    for p_kpa in P_vec_kPa:
        res_term = calcular_balanco_termico(F_kmol_h, xd, xb, zf, R_ref, q_ref, p_kpa * 1000.0, modelo, ctes_thermo, db_compostos, nome_v, nome_p)
        T_topo_C_P.append(res_term['T_topo_C'])
        T_fundo_C_P.append(res_term['T_fundo_C'])
        Qc_kW_P.append(res_term['Qc_kW'])

    q_vec = np.linspace(0.00, 1.00, 101)
    Qr_kW_q, n_estagios_q = [], []
    
    for q_val in q_vec:
        res_mc = calcular_mccabe_thiele_isobarico(xd, xb, zf, R_ref, q_val, murphree, tipo_murphree, P_total_Pa_atual, modelo, ctes_thermo, db_compostos, nome_v, nome_p)
        res_term = calcular_balanco_termico(F_kmol_h, xd, xb, zf, R_ref, q_val, P_total_Pa_atual, modelo, ctes_thermo, db_compostos, nome_v, nome_p)
        Qr_kW_q.append(res_term['Qr_kW'])
        n_estagios_q.append(res_mc['n_estagios'] if not res_mc['nao_convergiu'] else None)

    return {
        'R_vec': R_vec, 
        'n_estagios_R': n_estagios_R, 
        'Qc_kW_R': Qc_kW_R,
        'P_vec_kPa': P_vec_kPa, 
        'T_topo_C_P': T_topo_C_P, 
        'T_fundo_C_P': T_fundo_C_P, 
        'Qc_kW_P': Qc_kW_P,
        'q_vec': q_vec, 
        'Qr_kW_q': Qr_kW_q, 
        'n_estagios_q': n_estagios_q
    }

#---------------------------------------------- FORMATAÇÃO DE GRAFICOS PARA PDF ---------------------------------------#
def preparar_figura_para_pdf(fig, show_legend=True, height=450, width=720):
    
    fig_pdf = copy.deepcopy(fig)
    
    fig_pdf.update_layout(
        title=None,
        height=height,
        width=width,
        margin=dict(l=65, r=30, t=35, b=55),
        autosize=False
    )
    
    fig_pdf.update_xaxes(title_font=dict(size=12), tickfont=dict(size=10), title_standoff=12)
    fig_pdf.update_yaxes(title_font=dict(size=12), tickfont=dict(size=10), title_standoff=12)
    
    if show_legend:
        fig_pdf.update_layout(
            showlegend=True,
            legend=dict(
                orientation="h",
                yanchor="bottom",
                y=1.02,
                xanchor="center",
                x=0.5,
                font=dict(size=10)
            )
        )
    else:
        fig_pdf.update_layout(showlegend=False)
        
    return fig_pdf

#---------------------------------------------- GERACAO DO RELATORIO PDF -----------------------------------------------#
def gerar_pdf_relatorio(mistura, modelo, P_kPa, F, xd, zf, xb, R, q, murphree, tipo_murphree, 
                        mc_dados, balanco, fig_mc, fig_txy, fig_pxy, x1_grid, y1_isobarico, 
                        T_bolha_arr, Y_v_iso, Ptotal_iso, T_ref_K=300.0, sensib_dados=None, 
                        fig_sensib_R=None, fig_sensib_P=None, fig_sensib_q=None,
                        R_min=None, n_min=None):
  
    buffer = io.BytesIO()
    doc = SimpleDocTemplate(
        buffer, 
        pagesize=landscape(letter),
        rightMargin=36, leftMargin=36, topMargin=54, bottomMargin=36
    )
    story = []
    
    styles = getSampleStyleSheet()
    
    title_style = ParagraphStyle(
        'DocTitle', parent=styles['Heading1'],
        fontSize=18, leading=22, textColor=colors.HexColor('#003366'),
        alignment=1, spaceAfter=8
    )
    
    section_style = ParagraphStyle(
        'SectionHeader', parent=styles['Heading2'],
        fontSize=12, leading=16, textColor=colors.HexColor('#003366'),
        spaceBefore=10, spaceAfter=8
    )
    
    normal_style = ParagraphStyle('DocNormal', parent=styles['Normal'], fontSize=8.5, leading=11)
    author_style = ParagraphStyle('DocAuthor', parent=styles['Normal'], fontSize=8, leading=11, textColor=colors.HexColor('#333333'))

    
    story.append(Paragraph("<b>YuNno-Thermo | Relatório Técnico de Destilação Binária</b>", title_style))
    story.append(Spacer(1, 4))

    # ---------------------------------------------BLOCO DE AUTORIA E CONTATO -----------------------------------------#
    autores_data = [
        [
            Paragraph("<b>Desenvolvido por:</b>", author_style),
            Paragraph("<b>Brenno Fabricio Piva Ostroschi:</b> @brennopiva | brenno.piva@gmail.com", author_style),
            Paragraph("<b>Yuri Stankunavicius:</b> @yuri.stanks | yuri.ys.stank@gmail.com", author_style)
        ]
    ]
    t_autores = Table(autores_data, colWidths=[120, 340, 260])
    t_autores.setStyle(TableStyle([
        ('BACKGROUND', (0,0), (-1,-1), colors.HexColor('#EAECEE')),
        ('BOX', (0,0), (-1,-1), 0.5, colors.HexColor('#BDC3C7')),
        ('VALIGN', (0,0), (-1,-1), 'MIDDLE'),
        ('ALIGN', (0,0), (-1,-1), 'LEFT'),
        ('BOTTOMPADDING', (0,0), (-1,-1), 4),
        ('TOPPADDING', (0,0), (-1,-1), 4),
    ]))
    story.append(t_autores)
    story.append(Spacer(1, 10))

    # ---------------------------------------------------------SEÇÃO 1: PARÂMETROS ------------------------------------#
    story.append(Paragraph("<b>1. Parâmetros de Operação & Termodinâmica</b>", section_style))
    param_data = [
        [Paragraph(f"<b>Mistura Binária:</b> {mistura}", normal_style), Paragraph(f"<b>Vazão Alimentação (F):</b> {F:.1f} kmol/h", normal_style)],
        [Paragraph(f"<b>Modelo Termodinâmico:</b> {modelo}", normal_style), Paragraph(f"<b>Fração Alimentação (Zf):</b> {zf:.3f}", normal_style)],
        [Paragraph(f"<b>Pressão da Coluna:</b> {P_kPa:.2f} kPa", normal_style), Paragraph(f"<b>Fração Topo (Xd):</b> {xd:.3f}", normal_style)],
        [Paragraph(f"<b>Razão de Refluxo (R):</b> {R:.2f}", normal_style), Paragraph(f"<b>Fração Fundo (Xb):</b> {xb:.3f}", normal_style)],
        [Paragraph(f"<b>Estado Térmico (q):</b> {q:.2f}", normal_style), Paragraph(f"<b>Eficiência Murphree:</b> {murphree*100:.0f}% ({tipo_murphree})", normal_style)]
    ]
    t_param = Table(param_data, colWidths=[360, 360])
    t_param.setStyle(TableStyle([
        ('BACKGROUND', (0,0), (-1,-1), colors.HexColor('#F4F6F9')),
        ('BOX', (0,0), (-1,-1), 0.8, colors.HexColor('#CCCCCC')),
        ('VALIGN', (0,0), (-1,-1), 'MIDDLE'),
        ('BOTTOMPADDING', (0,0), (-1,-1), 3),
        ('TOPPADDING', (0,0), (-1,-1), 3),
    ]))
    story.append(t_param)
    story.append(Spacer(1, 8))

    # -------------------------------------------------------------SEÇÃO 2: MÉTRICAS -----------------------------------#    
    story.append(Paragraph("<b>2. Métricas de Dimensionamento e Cargas Térmicas</b>", section_style))
    r_min_str = f"{R_min:.4f}" if (R_min is not None and np.isfinite(R_min)) else "N/A"
    n_min_str = f"{n_min}" if n_min is not None else "N/A"
    est_alim = mc_dados.get('estagio_alimentacao')
    est_alim_str = f"{est_alim}º estágio" if est_alim is not None else "N/A"

    metric_data = [
        ["Estágios Reais (Alim.)", "Vazão Destilado (D)", "Temp. Topo", "Carga Condensador (Qc)", "Refluxo Mínimo (Rmín)"],
        [f"{mc_dados['n_estagios']} ({est_alim_str})", f"{balanco['D']:.2f} kmol/h", f"{balanco['T_topo_C']:.1f} °C", f"{balanco['Qc_kW']:.2f} kW", r_min_str],
        ["Eficiência Aplicada", "Vazão Resíduo (B)", "Temp. Fundo", "Carga Refervedor (Qr)", "Estágios Mínimos (Nmín)"],
        [f"{murphree*100:.0f}%", f"{balanco['B']:.2f} kmol/h", f"{balanco['T_fundo_C']:.1f} °C", f"{balanco['Qr_kW']:.2f} kW", n_min_str]
    ]
    t_metrics = Table(metric_data, colWidths=[144, 144, 144, 144, 144])
    t_metrics.setStyle(TableStyle([
        ('BACKGROUND', (0,0), (-1,0), colors.HexColor('#003366')),
        ('TEXTCOLOR', (0,0), (-1,0), colors.white),
        ('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'),
        ('BACKGROUND', (0,2), (-1,2), colors.HexColor('#003366')),
        ('TEXTCOLOR', (0,2), (-1,2), colors.white),
        ('FONTNAME', (0,2), (-1,2), 'Helvetica-Bold'),
        ('ALIGN', (0,0), (-1,-1), 'CENTER'),
        ('GRID', (0,0), (-1,-1), 0.5, colors.HexColor('#CCCCCC')),
        ('FONTSIZE', (0,0), (-1,-1), 8),
        ('BOTTOMPADDING', (0,0), (-1,-1), 4),
        ('TOPPADDING', (0,0), (-1,-1), 4),
    ]))
    story.append(t_metrics)

    def plotly_to_rl_image(fig, show_legend=True, width=720, height=450):
        fig_adj = preparar_figura_para_pdf(fig, show_legend=show_legend, height=height, width=width)
        img_bytes = fig_adj.to_image(format="png", width=width, height=height, scale=2)
        return RLImage(io.BytesIO(img_bytes), width=width, height=height)

    # -----------------------------------------------------3. McCABE-THIELE -------------------------------------------#
    story.append(PageBreak())
    story.append(Paragraph(f"<b>3. Diagrama de McCabe-Thiele Isobárico ({P_kPa:.1f} kPa)</b>", section_style))
    story.append(plotly_to_rl_image(fig_mc, show_legend=True))
    
    story.append(PageBreak())
    story.append(Paragraph("<b>3.1 Tabela de Dados Numéricos: Diagrama de McCabe-Thiele</b>", section_style))
    
    r_sl = R / (R + 1.0)
    b_sl = xd / (R + 1.0)
    x_q_pt = mc_dados['x_q'][1]
    m_lose = (mc_dados['y_q'][1] - xb) / (x_q_pt - xb) if abs(x_q_pt - xb) > 1e-6 else 1.0

    tb_mc_header = [["x1 (Líquido)", "y1 (Equilíbrio)", "y_Retificação", "y_Esgotamento"]]
    tb_mc_body = []
    for idx in range(len(x1_grid)):
        x_val = x1_grid[idx]
        y_eq = y1_isobarico[idx]
        y_ret = r_sl * x_val + b_sl if x_val >= x_q_pt else "-"
        y_esg = m_lose * (x_val - xb) + xb if x_val <= x_q_pt else "-"
        
        y_ret_str = f"{y_ret:.3f}" if isinstance(y_ret, float) else "-"
        y_esg_str = f"{y_esg:.3f}" if isinstance(y_esg, float) else "-"
        tb_mc_body.append([f"{x_val:.3f}", f"{y_eq:.3f}", y_ret_str, y_esg_str])

    t_mc_table = Table(tb_mc_header + tb_mc_body, colWidths=[180, 180, 180, 180])
    t_mc_table.setStyle(TableStyle([
        ('BACKGROUND', (0,0), (-1,0), colors.HexColor('#E60505')),
        ('TEXTCOLOR', (0,0), (-1,0), colors.white),
        ('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'),
        ('ALIGN', (0,0), (-1,-1), 'CENTER'),
        ('GRID', (0,0), (-1,-1), 0.5, colors.HexColor('#DDDDDD')),
        ('ROWBACKGROUNDS', (0,1), (-1,-1), [colors.white, colors.HexColor('#F9F9F9')]),
        ('FONTSIZE', (0,0), (-1,-1), 6),
        ('BOTTOMPADDING', (0,0), (-1,-1), 1.5),
        ('TOPPADDING', (0,0), (-1,-1), 1.5),
    ]))
    story.append(t_mc_table)

    # --------------------------------------------------------4. DIAGRAMA T-X-Y ---------------------------------------#
    story.append(PageBreak())
    story.append(Paragraph(f"<b>4. Diagrama T-x-y Isobárico ({P_kPa:.1f} kPa)</b>", section_style))
    story.append(plotly_to_rl_image(fig_txy, show_legend=True))
    
    story.append(PageBreak())
    story.append(Paragraph("<b>4.1 Tabela de Dados Numéricos: Diagrama T-x-y</b>", section_style))
    
    tb_txy_header = [["Fração Molar (x, y)", "Temperatura de Bolha (°C)", "Temperatura de Orvalho (°C)"]]
    tb_txy_body = []
    for idx in range(len(x1_grid)):
        x_val = x1_grid[idx]
        tb_c = T_bolha_arr[idx] - 273.15
        tb_txy_body.append([f"{x_val:.3f}", f"{tb_c:.2f}", f"{tb_c:.2f}"])

    t_txy_table = Table(tb_txy_header + tb_txy_body, colWidths=[240, 240, 240])
    t_txy_table.setStyle(TableStyle([
        ('BACKGROUND', (0,0), (-1,0), colors.HexColor('#0055FF')),
        ('TEXTCOLOR', (0,0), (-1,0), colors.white),
        ('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'),
        ('ALIGN', (0,0), (-1,-1), 'CENTER'),
        ('GRID', (0,0), (-1,-1), 0.5, colors.HexColor('#DDDDDD')),
        ('ROWBACKGROUNDS', (0,1), (-1,-1), [colors.white, colors.HexColor('#F9F9F9')]),
        ('FONTSIZE', (0,0), (-1,-1), 6),
        ('BOTTOMPADDING', (0,0), (-1,-1), 1.5),
        ('TOPPADDING', (0,0), (-1,-1), 1.5),
    ]))
    story.append(t_txy_table)

    # -------------------------------------------------- 5. DIAGRAMA P-X-Y ---------------------------------------------#
    story.append(PageBreak())
    story.append(Paragraph(f"<b>5. Diagrama P-x-y Isotérmico (Ref: {T_ref_K:.2f} K / {T_ref_K - 273.15:.2f} °C)</b>", section_style))
    story.append(plotly_to_rl_image(fig_pxy, show_legend=True))
    
    story.append(PageBreak())
    story.append(Paragraph("<b>5.1 Tabela de Dados Numéricos: Diagrama P-x-y</b>", section_style))
    
    tb_pxy_header = [["x1 (Líquido)", "y1 (Vapor)", "Pressão de Bolha Total (kPa)"]]
    tb_pxy_body = []
    for idx in range(len(x1_grid)):
        x_val = x1_grid[idx]
        y_val = Y_v_iso[idx]
        p_val = Ptotal_iso[idx] / 1000.0
        tb_pxy_body.append([f"{x_val:.3f}", f"{y_val:.3f}", f"{p_val:.2f}"])

    t_pxy_table = Table(tb_pxy_header + tb_pxy_body, colWidths=[240, 240, 240])
    t_pxy_table.setStyle(TableStyle([
        ('BACKGROUND', (0,0), (-1,0), colors.HexColor('#2E7D32')),
        ('TEXTCOLOR', (0,0), (-1,0), colors.white),
        ('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'),
        ('ALIGN', (0,0), (-1,-1), 'CENTER'),
        ('GRID', (0,0), (-1,-1), 0.5, colors.HexColor('#DDDDDD')),
        ('ROWBACKGROUNDS', (0,1), (-1,-1), [colors.white, colors.HexColor('#F9F9F9')]),
        ('FONTSIZE', (0,0), (-1,-1), 6),
        ('BOTTOMPADDING', (0,0), (-1,-1), 1.5),
        ('TOPPADDING', (0,0), (-1,-1), 1.5),
    ]))
    story.append(t_pxy_table)

    # ----------------------------------------- 6, 7 e 8. SENSIBILIDADES ---------------------------------------------#
    if sensib_dados is not None and fig_sensib_R is not None:
        
        story.append(PageBreak())
        story.append(Paragraph("<b>6. Análise de Sensibilidade: Razão de Refluxo (R)</b>", section_style))
        story.append(plotly_to_rl_image(fig_sensib_R, show_legend=True))
        
        story.append(PageBreak())
        story.append(Paragraph("<b>6.1 Tabela de Dados Numéricos: Sensibilidade R</b>", section_style))
        tb_r_header = [["Razão Refluxo (R)", "Nº Estágios Reais", "Carga Condensador Qc (kW)"]]
        tb_r_body = []
        for i in range(len(sensib_dados['R_vec'])):
            r_val = sensib_dados['R_vec'][i]
            n_est = sensib_dados['n_estagios_R'][i]
            qc = sensib_dados['Qc_kW_R'][i]
            tb_r_body.append([f"{r_val:.2f}", f"{n_est}" if n_est is not None else "Pinch", f"{qc:.2f}"])
            
        t_r_table = Table(tb_r_header + tb_r_body, colWidths=[240, 240, 240])
        t_r_table.setStyle(TableStyle([
            ('BACKGROUND', (0,0), (-1,0), colors.HexColor('#4169E1')),
            ('TEXTCOLOR', (0,0), (-1,0), colors.white),
            ('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'),
            ('ALIGN', (0,0), (-1,-1), 'CENTER'),
            ('GRID', (0,0), (-1,-1), 0.5, colors.HexColor('#DDDDDD')),
            ('ROWBACKGROUNDS', (0,1), (-1,-1), [colors.white, colors.HexColor('#F9F9F9')]),
            ('FONTSIZE', (0,0), (-1,-1), 6),
            ('BOTTOMPADDING', (0,0), (-1,-1), 1.5),
            ('TOPPADDING', (0,0), (-1,-1), 1.5),
        ]))
        story.append(t_r_table)

        
        story.append(PageBreak())
        story.append(Paragraph("<b>7. Análise de Sensibilidade: Pressão de Operação (P)</b>", section_style))
        story.append(plotly_to_rl_image(fig_sensib_P, show_legend=True))
        
        story.append(PageBreak())
        story.append(Paragraph("<b>7.1 Tabela de Dados Numéricos: Sensibilidade P</b>", section_style))
        tb_p_header = [["Pressão (kPa)", "Temp. Topo (°C)", "Temp. Fundo (°C)", "Carga Qc (kW)"]]
        tb_p_body = []
        for i in range(len(sensib_dados['P_vec_kPa'])):
            p_val = sensib_dados['P_vec_kPa'][i]
            tt = sensib_dados['T_topo_C_P'][i]
            tf = sensib_dados['T_fundo_C_P'][i]
            qc = sensib_dados['Qc_kW_P'][i]
            tb_p_body.append([f"{p_val:.1f}", f"{tt:.1f}", f"{tf:.1f}", f"{qc:.2f}"])
            
        t_p_table = Table(tb_p_header + tb_p_body, colWidths=[180, 180, 180, 180])
        t_p_table.setStyle(TableStyle([
            ('BACKGROUND', (0,0), (-1,0), colors.HexColor('#00BFFF')),
            ('TEXTCOLOR', (0,0), (-1,0), colors.white),
            ('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'),
            ('ALIGN', (0,0), (-1,-1), 'CENTER'),
            ('GRID', (0,0), (-1,-1), 0.5, colors.HexColor('#DDDDDD')),
            ('ROWBACKGROUNDS', (0,1), (-1,-1), [colors.white, colors.HexColor('#F9F9F9')]),
            ('FONTSIZE', (0,0), (-1,-1), 6),
            ('BOTTOMPADDING', (0,0), (-1,-1), 1.5),
            ('TOPPADDING', (0,0), (-1,-1), 1.5),
        ]))
        story.append(t_p_table)

        
        story.append(PageBreak())
        story.append(Paragraph("<b>8. Análise de Sensibilidade: Condição Térmica (q)</b>", section_style))
        story.append(plotly_to_rl_image(fig_sensib_q, show_legend=True))
        
        story.append(PageBreak())
        story.append(Paragraph("<b>8.1 Tabela de Dados Numéricos: Sensibilidade q</b>", section_style))
        tb_q_header = [["Parâmetro q", "Carga Refervedor Qr (kW)", "Nº Estágios Reais"]]
        tb_q_body = []
        for i in range(len(sensib_dados['q_vec'])):
            q_val = sensib_dados['q_vec'][i]
            qr = sensib_dados['Qr_kW_q'][i]
            n_est = sensib_dados['n_estagios_q'][i]
            tb_q_body.append([f"{q_val:.2f}", f"{qr:.2f}", f"{n_est}" if n_est is not None else "Pinch"])
            
        t_q_table = Table(tb_q_header + tb_q_body, colWidths=[240, 240, 240])
        t_q_table.setStyle(TableStyle([
            ('BACKGROUND', (0,0), (-1,0), colors.HexColor('#8B0000')),
            ('TEXTCOLOR', (0,0), (-1,0), colors.white),
            ('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'),
            ('ALIGN', (0,0), (-1,-1), 'CENTER'),
            ('GRID', (0,0), (-1,-1), 0.5, colors.HexColor('#DDDDDD')),
            ('ROWBACKGROUNDS', (0,1), (-1,-1), [colors.white, colors.HexColor('#F9F9F9')]),
            ('FONTSIZE', (0,0), (-1,-1), 6),
            ('BOTTOMPADDING', (0,0), (-1,-1), 1.5),
            ('TOPPADDING', (0,0), (-1,-1), 1.5),
        ]))
        story.append(t_q_table)

    doc.build(story)
    
    buffer.seek(0)
    return buffer.getvalue()

#---------------------------------------------- PACOTE ZIP (PDF + CSVs) ------------------------------------------------#
def gerar_pacote_zip(pdf_bytes, mistura, modelo, P_kPa, F, xd, zf, xb, R, q, murphree, tipo_murphree,
                     mc_dados, balanco, x1_grid, y1_isobarico, T_bolha_arr,
                     Y_v_iso, Ptotal_iso, T_ref_K=300.0, sensib_dados=None,
                     R_min=None, n_min=None):
    
    zip_buffer = io.BytesIO()
    sep = ";"

    def montar_csv(linhas):
        buf = io.StringIO()
        for linha in linhas:
            buf.write(sep.join(str(v) for v in linha) + "\n")
        return buf.getvalue().encode("utf-8-sig")

    
    linhas_metricas = [
        ["Parametro", "Valor", "Unidade"],
        ["Mistura Binaria", mistura, "-"],
        ["Modelo Termodinamico", modelo, "-"],
        ["Pressao da Coluna", f"{P_kPa:.2f}", "kPa"],
        ["Vazao Alimentacao (F)", f"{F:.2f}", "kmol/h"],
        ["Fracao Alimentacao (Zf)", f"{zf:.4f}", "-"],
        ["Fracao Topo (Xd)", f"{xd:.4f}", "-"],
        ["Fracao Fundo (Xb)", f"{xb:.4f}", "-"],
        ["Razao de Refluxo (R)", f"{R:.4f}", "-"],
        ["Razao de Refluxo Minima (Rmin)", f"{R_min:.4f}" if (R_min is not None and np.isfinite(R_min)) else "Pinch", "-"],
        ["Estado Termico (q)", f"{q:.4f}", "-"],
        ["Eficiencia Murphree", f"{murphree*100:.1f}", f"% ({tipo_murphree})"],
        ["Estagios Reais", mc_dados['n_estagios'], "estagios"],
        ["Estagios Minimos (Nmin)", n_min if n_min is not None else "N/A", "estagios"],
        ["Vazao Destilado (D)", f"{balanco['D']:.4f}", "kmol/h"],
        ["Vazao Residuo (B)", f"{balanco['B']:.4f}", "kmol/h"],
        ["Temperatura Topo", f"{balanco['T_topo_C']:.2f}", "°C"],
        ["Temperatura Fundo", f"{balanco['T_fundo_C']:.2f}", "°C"],
        ["Carga Condensador (Qc)", f"{balanco['Qc_kW']:.4f}", "kW"],
        ["Carga Refervedor (Qr)", f"{balanco['Qr_kW']:.4f}", "kW"],
    ]

    r_sl = R / (R + 1.0)
    b_sl = xd / (R + 1.0)
    x_q_pt = mc_dados['x_q'][1]
    m_lose = (mc_dados['y_q'][1] - xb) / (x_q_pt - xb) if abs(x_q_pt - xb) > 1e-6 else 1.0

    linhas_txy_mc = [["x1 (Liquido)", "y1 (Vapor Equilibrio)", "T_Equilibrio (°C)", "y_Retificacao", "y_Esgotamento"]]
    for idx in range(len(x1_grid)):
        x_val = x1_grid[idx]
        y_eq = y1_isobarico[idx]
        tb_c = T_bolha_arr[idx] - 273.15
        y_ret = f"{(r_sl * x_val + b_sl):.4f}" if x_val >= x_q_pt else ""
        y_esg = f"{(m_lose * (x_val - xb) + xb):.4f}" if x_val <= x_q_pt else ""
        linhas_txy_mc.append([f"{x_val:.4f}", f"{y_eq:.4f}", f"{tb_c:.2f}", y_ret, y_esg])

   
    linhas_pxy = [["x1 (Liquido)", "y1 (Vapor)", f"Pressao Total a {T_ref_K:.2f} K (kPa)"]]
    for idx in range(len(x1_grid)):
        linhas_pxy.append([f"{x1_grid[idx]:.4f}", f"{Y_v_iso[idx]:.4f}", f"{Ptotal_iso[idx] / 1000.0:.4f}"])

    with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zf_out:
        zf_out.writestr(f"Relatorio_Destilacao_{mistura}_{P_kPa:.0f}kPa.pdf", pdf_bytes)
        zf_out.writestr("01_Parametros_e_Metricas.csv", montar_csv(linhas_metricas))
        zf_out.writestr("02_Equilibrio_TXY_e_McCabe_Thiele.csv", montar_csv(linhas_txy_mc))
        zf_out.writestr("03_Equilibrio_PXY_Isotermico.csv", montar_csv(linhas_pxy))


        if sensib_dados is not None:
            linhas_sens_R = [["Razao_Refluxo_R", "Estagios_Reais", "Carga_Condensador_Qc_kW"]]
            for i in range(len(sensib_dados['R_vec'])):
                est_r = sensib_dados['n_estagios_R'][i]
                linhas_sens_R.append([
                    f"{sensib_dados['R_vec'][i]:.4f}",
                    est_r if est_r is not None else "Pinch",
                    f"{sensib_dados['Qc_kW_R'][i]:.4f}"
                ])
            zf_out.writestr("04_Sensibilidade_Refluxo_R.csv", montar_csv(linhas_sens_R))

            linhas_sens_P = [["Pressao_kPa", "Temp_Topo_C", "Temp_Fundo_C", "Carga_Condensador_Qc_kW"]]
            for i in range(len(sensib_dados['P_vec_kPa'])):
                linhas_sens_P.append([
                    f"{sensib_dados['P_vec_kPa'][i]:.2f}",
                    f"{sensib_dados['T_topo_C_P'][i]:.2f}",
                    f"{sensib_dados['T_fundo_C_P'][i]:.2f}",
                    f"{sensib_dados['Qc_kW_P'][i]:.4f}"
                ])
            zf_out.writestr("05_Sensibilidade_Pressao_P.csv", montar_csv(linhas_sens_P))

            linhas_sens_q = [["Parametro_q", "Carga_Refervedor_Qr_kW", "Estagios_Reais"]]
            for i in range(len(sensib_dados['q_vec'])):
                est_q = sensib_dados['n_estagios_q'][i]
                linhas_sens_q.append([
                    f"{sensib_dados['q_vec'][i]:.4f}",
                    f"{sensib_dados['Qr_kW_q'][i]:.4f}",
                    est_q if est_q is not None else "Pinch"
                ])
            zf_out.writestr("06_Sensibilidade_Alimentacao_q.csv", montar_csv(linhas_sens_q))

    zip_buffer.seek(0)
    return zip_buffer.getvalue()