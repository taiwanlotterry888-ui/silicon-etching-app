# -*- coding: utf-8 -*-
"""
Silicon Etching Process Evaluator (v2 - Physics-based)
--------------------------------------------------------
矽（Silicon）蝕刻製程評估系統 - 動力學模型版

本版本核心改進：濕蝕刻與乾蝕刻速率不再只是查表，而是採用
Arrhenius 動力學方程式，依據使用者輸入的「溫度」實際計算蝕刻速率。

資料來源總覽（詳細出處見程式內各段落註解，以及頁面⑦「參考文獻」）：
 1. 課堂 PDF《Etching》(Kovacs, Maluf & Peterson 1998; Williams & Muller)
    -- 材料選擇比查詢表格
 2. Seidel, H. et al. (1990). J. Electrochem. Soc. 137, 3612 / 3626
    -- KOH 動力學參數 (Ea, R0)，硼掺杂 etch-stop
 3. Shikida, M. et al. (2000). Sens. Actuators A 80, 179-188
    -- KOH / TMAH 動力學參數
 4. Tabata, O. et al. -- TMAH 動力學參數
 5. Dutta, S. et al. (2011). Microsystem Technologies
    -- KOH/TMAH/EDP 活化能比較 (Si(110))
 6. Gosálvez, M.A., Zubel, I., Viinikka, E. (2015). "Wet Etching of
    Silicon", in Handbook of Silicon Based MEMS Materials and
    Technologies (2nd ed.), Ch.22, Table 22.3, Eq.(22.3)
    -- 彙整上述多篇 KOH/TMAH 動力學數據的二次文獻
 7. Ohring, M. (2002). Materials Science of Thin Films (2nd ed.),
    Ch.5, Eq.(5-9a)(5-9b)
    -- 乾蝕刻（電漿氟原子化學反應）動力學模型

限制與待查事項（誠實列出，供之後查證）：
 - KOH/TMAH 的「濃度」相依性目前沒有做成公式，因為 Seidel 浓度公式
   換算 wt%→mol/L 所需的密度多項式系數 (a0,a1,a2) 追不到原始出處，
   因此本版本 KOH/TMAH/EDP 均為「固定濃度、僅溫度可調」的模型。
 - EDP 的 R0 前置因子文獻未給出數值，是用課堂 PDF 的單一實測點
   (115°C, 0.75 µm/min) 反推得出，非文獻原始數值，UI 會特別註明。

部署方式：
1. 把整個資料夾 push 到 GitHub
2. 到 https://share.streamlit.io 用該 GitHub repo 部署
3. 主檔案填 app.py，即可產生網址
"""

import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

st.set_page_config(page_title="Silicon Etching Process Evaluator",
                    layout="wide")

KB_EV = 8.617333262e-5  # Boltzmann constant, eV/K

# =========================================================
# 1. 濕蝕刻動力學資料庫 (Arrhenius: R = R0 * exp(-Ea/kB T))
#    單位：R0 以 µm/h 儲存（多數文獻的原始單位），計算後換算成 µm/min
# =========================================================
wet_kinetics = pd.DataFrame([
    # ---- KOH ----
    {"Etchant": "KOH", "Concentration": "34 wt%", "Plane": "{100}",
     "Ea_eV": 0.61, "R0_um_h": 3.10e10, "Source": "Seidel et al. 1990"},
    {"Etchant": "KOH", "Concentration": "34 wt%", "Plane": "{110}",
     "Ea_eV": 0.60, "R0_um_h": 3.66e10, "Source": "Seidel et al. 1990"},
    {"Etchant": "KOH", "Concentration": "34 wt%", "Plane": "{100}",
     "Ea_eV": 0.62, "R0_um_h": 4.44e10, "Source": "Shikida et al. 2000"},
    {"Etchant": "KOH", "Concentration": "34 wt%", "Plane": "{110}",
     "Ea_eV": 0.60, "R0_um_h": 5.17e10, "Source": "Shikida et al. 2000"},
    {"Etchant": "KOH", "Concentration": "34 wt%", "Plane": "{111}",
     "Ea_eV": 0.48, "R0_um_h": 7.31e6, "Source": "Shikida et al. 2000"},
    {"Etchant": "KOH", "Concentration": "35 wt%", "Plane": "{100}",
     "Ea_eV": 0.62, "R0_um_h": 4.16e10, "Source": "Tan et al."},
    {"Etchant": "KOH", "Concentration": "35 wt%", "Plane": "{111}",
     "Ea_eV": 0.52, "R0_um_h": 1.84e7, "Source": "Tan et al."},
    {"Etchant": "KOH", "Concentration": "37 wt%", "Plane": "{100}",
     "Ea_eV": 0.56, "R0_um_h": 6.06e9, "Source": "Wind et al."},
    {"Etchant": "KOH", "Concentration": "37 wt%", "Plane": "{110}",
     "Ea_eV": 0.56, "R0_um_h": 1.21e10, "Source": "Wind et al."},
    {"Etchant": "KOH", "Concentration": "37 wt%", "Plane": "{111}",
     "Ea_eV": 0.55, "R0_um_h": 3.93e8, "Source": "Wind et al."},
    # ---- TMAH ----
    {"Etchant": "TMAH", "Concentration": "22 wt%", "Plane": "{100}",
     "Ea_eV": 0.65, "R0_um_h": 6.40e10, "Source": "Tabata et al."},
    {"Etchant": "TMAH", "Concentration": "22 wt%", "Plane": "{110}",
     "Ea_eV": 0.39, "R0_um_h": 1.83e7, "Source": "Tabata et al."},
    {"Etchant": "TMAH", "Concentration": "22 wt%", "Plane": "{111}",
     "Ea_eV": 0.71, "R0_um_h": 1.56e10, "Source": "Tabata et al."},
    {"Etchant": "TMAH", "Concentration": "25 wt%", "Plane": "{100}",
     "Ea_eV": 0.65, "R0_um_h": 5.51e10, "Source": "Shikida et al."},
    {"Etchant": "TMAH", "Concentration": "25 wt%", "Plane": "{110}",
     "Ea_eV": 0.64, "R0_um_h": 7.05e10, "Source": "Shikida et al."},
    {"Etchant": "TMAH", "Concentration": "25 wt%", "Plane": "{111}",
     "Ea_eV": 0.73, "R0_um_h": 2.55e10, "Source": "Shikida et al."},
])

# ---- EDP：文獻只給 Ea，R0 用課堂 PDF 單點資料反推 ----
EDP_EA = 0.34  # eV, Dutta et al. 2011 (Si(110)), 無特定濃度標示
EDP_REF_T_C = 115.0        # PDF General Characteristics 表：EDP 溫度
EDP_REF_RATE_UM_MIN = 0.75  # PDF 表：EDP 蝕刻速率 (µm/min)


def edp_r0_um_h():
    """由課堂 PDF 單一資料點反推 EDP 的 Arrhenius 前置因子 R0 (µm/h)。"""
    T_K = EDP_REF_T_C + 273.15
    R_um_h = EDP_REF_RATE_UM_MIN * 60.0
    return R_um_h / np.exp(-EDP_EA / (KB_EV * T_K))


def arrhenius_rate_um_min(Ea_eV, R0_um_h, T_C):
    """Arrhenius: R = R0 * exp(-Ea/kB T)，回傳 µm/min。"""
    T_K = T_C + 273.15
    R_um_h = R0_um_h * np.exp(-Ea_eV / (KB_EV * T_K))
    return R_um_h / 60.0


# =========================================================
# 2. 乾蝕刻（電漿）動力學模型
#    出處：Ohring, M. (2002) Materials Science of Thin Films,
#    2nd ed., Ch.5, Eq.(5-9a)(5-9b)
#    Re = A * sqrt(T) * C_F * exp(-Ea/kB T)   單位：Å/min
#    T: 基板溫度 (K)；C_F: 氟原子濃度 (atoms/cm^3)
# =========================================================
DRY_PARAMS = {
    "Si":   {"A": 2.91e-12, "Ea_eV": 0.108},
    "SiO2": {"A": 6.14e-13, "Ea_eV": 0.163},
}


def dry_etch_rate_A_min(material, T_C, C_F):
    T_K = T_C + 273.15
    p = DRY_PARAMS[material]
    return p["A"] * np.sqrt(T_K) * C_F * np.exp(-p["Ea_eV"] / (KB_EV * T_K))


# =========================================================
# 3. 蝕刻劑 × 材料 選擇比查詢資料庫（沿用課堂 PDF 表格）
# =========================================================
etchant_material_rate = pd.DataFrame([
    {"Etchant": "Concentrated HF (49%)", "Type": "Wet",
     "TargetMaterial": "Silicon oxides",
     "Polysilicon_n+": 0, "Polysilicon_undoped": None,
     "SiO2": 2300, "SiNitride": 14, "PSG": 3600,
     "Aluminum": 4.2, "Titanium": 1000, "Photoresist": 0},
    {"Etchant": "25:1 HF:H2O", "Type": "Wet",
     "TargetMaterial": "Silicon oxides",
     "Polysilicon_n+": 0, "Polysilicon_undoped": 0,
     "SiO2": 9.7, "SiNitride": 0.6, "PSG": 150,
     "Aluminum": None, "Titanium": None, "Photoresist": 0},
    {"Etchant": "5:1 BHF", "Type": "Wet",
     "TargetMaterial": "Silicon oxides",
     "Polysilicon_n+": 9, "Polysilicon_undoped": 2,
     "SiO2": 100, "SiNitride": 0.9, "PSG": 440,
     "Aluminum": 140, "Titanium": 1000, "Photoresist": 0},
    {"Etchant": "Silicon etchant (126HNO3:60H2O:5NH4F)", "Type": "Wet",
     "TargetMaterial": "Silicon",
     "Polysilicon_n+": 310, "Polysilicon_undoped": 100,
     "SiO2": 9, "SiNitride": 0.2, "PSG": 170,
     "Aluminum": 400, "Titanium": 300, "Photoresist": 0},
    {"Etchant": "Aluminum etchant (16H3PO4:1HNO3:2H2)", "Type": "Wet",
     "TargetMaterial": "Aluminum",
     "Polysilicon_n+": 1, "Polysilicon_undoped": None,
     "SiO2": 12, "SiNitride": 0, "PSG": 1,
     "Aluminum": 660, "Titanium": 0, "Photoresist": 0},
    {"Etchant": "Titanium etchant (20H2O:1H2O2:1HF)", "Type": "Wet",
     "TargetMaterial": "Titanium/Metals/organics",
     "Polysilicon_n+": 1.2, "Polysilicon_undoped": None,
     "SiO2": 0, "SiNitride": 0.8, "PSG": 210,
     "Aluminum": 10, "Titanium": 880, "Photoresist": 0},
    {"Etchant": "Piranha (50H2SO4:1H2O2)", "Type": "Wet",
     "TargetMaterial": "Metals/organics (cleaning)",
     "Polysilicon_n+": 0, "Polysilicon_undoped": 0,
     "SiO2": 0, "SiNitride": 0, "PSG": 0,
     "Aluminum": 180, "Titanium": 240, "Photoresist": 10},
    {"Etchant": "Acetone (CH3COOH)", "Type": "Wet",
     "TargetMaterial": "Photoresist",
     "Polysilicon_n+": 0, "Polysilicon_undoped": 0,
     "SiO2": 0, "SiNitride": 0, "PSG": 0,
     "Aluminum": 0, "Titanium": 0, "Photoresist": 4000},
    {"Etchant": "CF4+CHF3+He, 450W", "Type": "Dry",
     "TargetMaterial": "Silicon oxides/nitrides",
     "Polysilicon_n+": 190, "Polysilicon_undoped": 210,
     "SiO2": 470, "SiNitride": 180, "PSG": 620,
     "Aluminum": None, "Titanium": 1000, "Photoresist": 220},
    {"Etchant": "SF6+He, 100W", "Type": "Dry",
     "TargetMaterial": "Thin silicon nitrides",
     "Polysilicon_n+": 73, "Polysilicon_undoped": 67,
     "SiO2": 31, "SiNitride": 32, "PSG": 61,
     "Aluminum": None, "Titanium": 1000, "Photoresist": 69},
    {"Etchant": "SF6, 12.5W", "Type": "Dry",
     "TargetMaterial": "Thin silicon nitrides",
     "Polysilicon_n+": 170, "Polysilicon_undoped": 280,
     "SiO2": 110, "SiNitride": 280, "PSG": 140,
     "Aluminum": None, "Titanium": 1000, "Photoresist": 310},
    {"Etchant": "O2, 400W", "Type": "Dry",
     "TargetMaterial": "Ashing photoresist",
     "Polysilicon_n+": 0, "Polysilicon_undoped": 0,
     "SiO2": 0, "SiNitride": 0, "PSG": 0,
     "Aluminum": 0, "Titanium": 0, "Photoresist": 340},
])

MATERIAL_COLS = ["Polysilicon_n+", "Polysilicon_undoped", "SiO2",
                  "SiNitride", "PSG", "Aluminum", "Titanium", "Photoresist"]

ANGLE_111 = 54.7  # {111} 面角度 (度)，異向性蝕刻幾何常數

# =========================================================
# Sidebar 導覽
# =========================================================
st.sidebar.title("🔬 Silicon Etching Evaluator")
page = st.sidebar.radio(
    "選擇功能模組",
    [
        "① 首頁說明",
        "② 濕蝕刻動力學計算器（KOH / TMAH / EDP）",
        "③ 乾蝕刻化學動力學計算器（F原子模型）",
        "④ 蝕刻劑材料選擇比查詢",
        "⑤ 異向性蝕刻剖面模擬器",
        "⑥ Mask 材料建議",
        "⑦ 參考文獻與資料來源",
    ],
)

# =========================================================
# ① 首頁
# =========================================================
if page == "① 首頁說明":
    st.title("矽蝕刻製程評估系統（動力學模型版）")
    st.markdown("""
本版本的核心改進：**濕蝕刻與乾蝕刻速率不再是單純查表**，而是採用
Arrhenius 動力學方程式 $R = R_0 \\cdot e^{-E_a/k_BT}$，
依照使用者輸入的**實際溫度**去計算蝕刻速率，所有係數（$E_a$、$R_0$）
皆標註文獻出處，可在頁面⑦「參考文獻」逐條查核。

### 模組總覽
| 模組 | 物理模型 | 可調參數 |
|---|---|---|
| 濕蝕刻動力學計算器 | Arrhenius（固定濃度） | 蝕刻劑、晶面、文獻來源、溫度 |
| 乾蝕刻化學動力學計算器 | Arrhenius + 氟原子濃度模型 | 材料、溫度、氟原子濃度 |
| 蝕刻劑材料選擇比查詢 | 查表（課堂 PDF 原始資料） | 蝕刻劑、材料 |
| 異向性蝕刻剖面模擬器 | 幾何模型（{111}面 54.7°） | 開口寬度、深度 |
| Mask 材料建議 | 查表 + 選擇比計算 | 蝕刻劑 |

### 已知限制（誠實聲明）
- KOH / TMAH 的「濃度」目前**不是**可連續調整的變數，僅能在文獻報告過的
  固定濃度中選擇（例如 KOH 34/35/37 wt%），因為 wt%→mol/L 換算所需的
  密度公式係數查無可靠原始出處，強行套用會誤導使用者。
- EDP 的前置因子 $R_0$ 是由課堂 PDF 單點資料反推，非文獻直接數值。
""")

# =========================================================
# ② 濕蝕刻動力學計算器
# =========================================================
elif page == "② 濕蝕刻動力學計算器（KOH / TMAH / EDP）":
    st.header("濕蝕刻動力學計算器")
    st.latex(r"R = R_0 \cdot e^{-E_a / k_B T}")
    st.caption("固定濃度模型：濃度已內建於所選文獻資料組中，僅溫度可調整。")

    etchant_choice = st.selectbox("選擇蝕刻劑", ["KOH", "TMAH", "EDP"])

    if etchant_choice in ["KOH", "TMAH"]:
        sub = wet_kinetics[wet_kinetics["Etchant"] == etchant_choice]
        # 讓使用者選擇 濃度+晶面+文獻來源 的組合
        options = sub.apply(
            lambda r: f"{r['Concentration']} | {r['Plane']} | {r['Source']}",
            axis=1).tolist()
        choice = st.selectbox("選擇文獻資料組（濃度 / 晶面 / 研究團隊）", options)
        idx = options.index(choice)
        row = sub.iloc[idx]

        col1, col2 = st.columns(2)
        with col1:
            st.metric("活化能 Ea (eV)", f"{row['Ea_eV']}")
            st.metric("濃度", row["Concentration"])
        with col2:
            st.metric("前置因子 R0 (µm/h)", f"{row['R0_um_h']:.3e}")
            st.metric("晶面", row["Plane"])
        st.caption(f"資料來源：{row['Source']}（詳見頁面⑦）")

        T_C = st.slider("蝕刻溫度 (°C)", 20.0, 130.0, 80.0, 1.0)
        rate = arrhenius_rate_um_min(row["Ea_eV"], row["R0_um_h"], T_C)
        st.success(f"預估蝕刻速率：**{rate:.3f} µm/min**"
                    f"（{rate*60:.2f} µm/h）　＠ {T_C:.0f}°C")

        # Arrhenius 曲線圖
        T_range = np.linspace(20, 130, 200)
        R_range = [arrhenius_rate_um_min(row["Ea_eV"], row["R0_um_h"], t)
                   for t in T_range]
        fig, ax = plt.subplots(figsize=(6, 3.5))
        ax.plot(T_range, R_range, 'b-')
        ax.axvline(T_C, color='red', linestyle='--', alpha=0.6)
        ax.plot([T_C], [rate], 'ro')
        ax.set_xlabel("Temperature (°C)")
        ax.set_ylabel("Etch rate (µm/min)")
        ax.set_yscale('log')
        ax.set_title(f"{etchant_choice} Arrhenius plot ({row['Plane']}, {row['Concentration']})")
        ax.grid(alpha=0.3)
        st.pyplot(fig)

        st.markdown("---")
        st.subheader("目標深度 → 預估時間")
        target_depth_um = st.number_input("目標蝕刻深度 (µm)", min_value=0.0,
                                           value=100.0, step=10.0)
        if rate > 0:
            st.info(f"預估所需時間：約 **{target_depth_um/rate:.1f} 分鐘**")

    else:  # EDP
        st.warning("⚠️ EDP 的前置因子 R0 並非文獻直接給出，"
                    "而是用課堂 PDF 的單點資料"
                    f"（{EDP_REF_T_C:.0f}°C, {EDP_REF_RATE_UM_MIN} µm/min）"
                    "反推計算，準確度僅供參考。")
        R0 = edp_r0_um_h()
        col1, col2 = st.columns(2)
        with col1:
            st.metric("活化能 Ea (eV)", f"{EDP_EA}")
            st.caption("來源：Dutta et al. 2011（Si(110)）")
        with col2:
            st.metric("反推前置因子 R0 (µm/h)", f"{R0:.3e}")
            st.caption("由課堂 PDF 單點資料反推")

        T_C = st.slider("蝕刻溫度 (°C)", 60.0, 130.0, 115.0, 1.0)
        rate = arrhenius_rate_um_min(EDP_EA, R0, T_C)
        st.success(f"預估蝕刻速率：**{rate:.3f} µm/min**　＠ {T_C:.0f}°C")

        T_range = np.linspace(60, 130, 200)
        R_range = [arrhenius_rate_um_min(EDP_EA, R0, t) for t in T_range]
        fig, ax = plt.subplots(figsize=(6, 3.5))
        ax.plot(T_range, R_range, 'g-')
        ax.axvline(T_C, color='red', linestyle='--', alpha=0.6)
        ax.plot([T_C], [rate], 'ro')
        ax.plot([EDP_REF_T_C], [EDP_REF_RATE_UM_MIN], 'k*', markersize=12,
                label='PDF 實測點')
        ax.set_xlabel("Temperature (°C)")
        ax.set_ylabel("Etch rate (µm/min)")
        ax.legend()
        ax.set_title("EDP Arrhenius plot（R0反推版）")
        ax.grid(alpha=0.3)
        st.pyplot(fig)

    st.markdown("---")
    st.subheader("完整濕蝕刻動力學資料表（KOH / TMAH）")
    st.dataframe(wet_kinetics, use_container_width=True)

# =========================================================
# ③ 乾蝕刻化學動力學計算器
# =========================================================
elif page == "③ 乾蝕刻化學動力學計算器（F原子模型）":
    st.header("乾蝕刻化學動力學計算器")
    st.latex(r"R_e = A \cdot T^{1/2} \cdot C_F \cdot e^{-E_a/k_BT}")
    st.caption("出處：Ohring, M. (2002) Materials Science of Thin Films, "
               "2nd ed., Eq.(5-9a)(5-9b)。R_e 單位 Å/min，T 為基板溫度(K)，"
               "C_F 為氟原子濃度 (atoms/cm³)。")

    col1, col2 = st.columns(2)
    with col1:
        T_C = st.slider("基板溫度 (°C)", 0.0, 200.0, 25.0, 5.0)
    with col2:
        C_F_exp = st.slider("氟原子濃度 C_F 指數 (10^x atoms/cm³)",
                             13.0, 17.0, 15.0, 0.1)
        C_F = 10 ** C_F_exp

    rate_si = dry_etch_rate_A_min("Si", T_C, C_F)
    rate_sio2 = dry_etch_rate_A_min("SiO2", T_C, C_F)

    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Si 蝕刻速率", f"{rate_si:.2f} Å/min",
                   help=f"= {rate_si/10:.3f} nm/min")
    with col2:
        st.metric("SiO2 蝕刻速率", f"{rate_sio2:.4f} Å/min",
                   help=f"= {rate_sio2/10:.4f} nm/min")
    with col3:
        if rate_sio2 > 0:
            st.metric("Si : SiO2 選擇比", f"{rate_si/rate_sio2:.1f} : 1")
        else:
            st.metric("Si : SiO2 選擇比", "∞（SiO2幾乎不蝕刻）")

    st.markdown("---")
    st.subheader("溫度 vs. 選擇比 關係圖")
    st.caption("Si 的活化能 (0.108 eV) 低於 SiO2 (0.163 eV)，"
               "代表溫度升高時 SiO2 速率上升更快，選擇比會隨溫度升高而下降。")
    T_range = np.linspace(0, 200, 100)
    sel_range = [dry_etch_rate_A_min("Si", t, C_F) /
                 max(dry_etch_rate_A_min("SiO2", t, C_F), 1e-12)
                 for t in T_range]
    fig, ax = plt.subplots(figsize=(6, 3.5))
    ax.plot(T_range, sel_range, 'purple')
    ax.axvline(T_C, color='red', linestyle='--', alpha=0.6)
    ax.set_xlabel("Substrate Temperature (°C)")
    ax.set_ylabel("Si : SiO2 Selectivity")
    ax.set_title(f"Selectivity vs. Temperature (C_F = 10^{C_F_exp:.1f} atoms/cm³)")
    ax.grid(alpha=0.3)
    st.pyplot(fig)

    st.markdown("---")
    st.subheader("目標深度 → 預估時間")
    target_nm = st.number_input("目標蝕刻深度 (nm)", min_value=0.0,
                                 value=500.0, step=50.0)
    rate_nm_min = rate_si / 10.0
    if rate_nm_min > 0:
        st.info(f"以 Si 蝕刻速率估算，預估所需時間：約 "
                f"**{target_nm/rate_nm_min:.1f} 分鐘**")

# =========================================================
# ④ 蝕刻劑材料選擇比查詢
# =========================================================
elif page == "④ 蝕刻劑材料選擇比查詢":
    st.header("蝕刻劑 × 材料 蝕刻速率／選擇比查詢")
    st.caption("資料來源：課堂 PDF《Etching》- "
               "Comparison of Etch Rates for Selected Etchants and Target Materials")

    tab1, tab2 = st.tabs(["依蝕刻劑查詢", "依材料反查"])

    with tab1:
        etch_choice = st.selectbox("選擇蝕刻劑",
                                    etchant_material_rate["Etchant"].tolist())
        row = etchant_material_rate[
            etchant_material_rate["Etchant"] == etch_choice].iloc[0]
        st.write(f"**目標材料（原設計用途）：** {row['TargetMaterial']}　"
                 f"**類型：** {row['Type']}")

        rates = {m: row[m] for m in MATERIAL_COLS}
        rates_df = pd.DataFrame(
            [{"材料": k, "蝕刻速率 (nm/min)": ("無資料" if pd.isna(v) else v)}
             for k, v in rates.items()])
        st.dataframe(rates_df, use_container_width=True)

        st.markdown("#### 計算兩種材料間的選擇比")
        m1 = st.selectbox("材料 A（要保留的，如 Mask）", MATERIAL_COLS, index=2)
        m2 = st.selectbox("材料 B（要蝕刻的）", MATERIAL_COLS, index=0)
        v1, v2 = row[m1], row[m2]
        if pd.isna(v1) or pd.isna(v2):
            st.warning("此蝕刻劑對所選材料之一無資料，無法計算選擇比。")
        elif v1 == 0:
            st.success(f"材料 A（{m1}）幾乎不被蝕刻 → 選擇比視為極高（理想 Mask 材料）")
        else:
            st.info(f"選擇比 (B:A) ≈ **{v2/v1:.1f} : 1**"
                    f"　（B={m2} 速率 {v2} nm/min ／ A={m1} 速率 {v1} nm/min）")

    with tab2:
        material_choice = st.selectbox("選擇要蝕刻的材料", MATERIAL_COLS)
        sub2 = etchant_material_rate[["Etchant", "Type", material_choice]].copy()
        sub2 = sub2.rename(columns={material_choice: "EtchRate_nm_min"})
        sub2 = sub2.dropna(subset=["EtchRate_nm_min"]).sort_values(
            "EtchRate_nm_min", ascending=False)
        st.write(f"針对 **{material_choice}**，各蝕刻劑速率（由快到慢）：")
        st.dataframe(sub2, use_container_width=True)
        if len(sub2) > 0:
            st.success(f"蝕刻速度最快：**{sub2.iloc[0]['Etchant']}** "
                        f"({sub2.iloc[0]['EtchRate_nm_min']} nm/min)")
            st.info(f"蝕刻速度最慢（最適合當作對此材料無害的 Mask）："
                    f"**{sub2.iloc[-1]['Etchant']}** "
                    f"({sub2.iloc[-1]['EtchRate_nm_min']} nm/min)")

    st.markdown("---")
    st.subheader("完整資料表")
    st.dataframe(etchant_material_rate, use_container_width=True)

# =========================================================
# ⑤ 異向性蝕刻剖面模擬器
# =========================================================
elif page == "⑤ 異向性蝕刻剖面模擬器":
    st.header("異向性蝕刻（Anisotropic Etching）剖面模擬")
    st.caption("模擬 <100> 矽晶片以 KOH/TMAH/EDP 進行 orientation-dependent "
               "etching (ODE) 的剖面。{111}面角度 54.7° 為晶格幾何常數。")

    col1, col2, col3 = st.columns(3)
    with col1:
        opening_um = st.number_input("Mask 開口寬度 (µm)", min_value=1.0,
                                      value=100.0, step=5.0)
    with col2:
        depth_um = st.number_input("蝕刻深度 (µm)", min_value=1.0,
                                    value=50.0, step=5.0)
    with col3:
        wafer_thickness_um = st.number_input("晶片厚度 (µm，選填)",
                                              min_value=0.0, value=0.0, step=10.0)

    mode = st.radio("蝕刻模式", ["異向性 (ODE, 54.7° {111}面)", "等向性 (Isotropic)"])

    if mode.startswith("異向性"):
        shrink_per_side = depth_um / np.tan(np.radians(ANGLE_111))
        bottom_width = opening_um - 2 * shrink_per_side
        if bottom_width <= 0:
            st.error("⚠️ 蝕刻深度過深，已形成 V 形溝槏（無平坦底部）！"
                      f" 理論最大深度 ≈ {opening_um/2*np.tan(np.radians(ANGLE_111)):.1f} µm")
            bottom_width = 0
            v_groove = True
        else:
            v_groove = False
            st.success(f"底部寬度 ≈ **{bottom_width:.1f} µm**"
                        f"（每側內縮 {shrink_per_side:.1f} µm）")
    else:
        shrink_per_side = depth_um
        bottom_width = opening_um
        v_groove = False
        st.info(f"等向性蝕刻：Mask 下方 undercut ≈ **{shrink_per_side:.1f} µm**")

    if wafer_thickness_um > 0 and depth_um >= wafer_thickness_um:
        st.warning("⚠️ 蝕刻深度已達到或超過晶片厚度，可能已貫穿晶片！")

    fig, ax = plt.subplots(figsize=(6, 4))
    top_y, bottom_y = 0, -depth_um
    half_open = opening_um / 2

    if mode.startswith("異向性"):
        if v_groove:
            apex_depth = half_open * np.tan(np.radians(ANGLE_111))
            xs, ys = [-half_open, 0, half_open], [top_y, -apex_depth, top_y]
        else:
            half_bottom = bottom_width / 2
            xs = [-half_open, -half_bottom, half_bottom, half_open]
            ys = [top_y, bottom_y, bottom_y, top_y]
        ax.plot(xs, ys, 'b-', linewidth=2)
        ax.fill(xs, ys, alpha=0.3, color='deepskyblue')
    else:
        half_bottom = half_open + shrink_per_side
        xs = [-half_open, -half_bottom, half_bottom, half_open]
        ys = [top_y, bottom_y, bottom_y, top_y]
        ax.plot(xs, ys, 'g-', linewidth=2)
        ax.fill(xs, ys, alpha=0.3, color='lightgreen')

    mask_w = opening_um * 1.5
    ax.plot([-mask_w/2, -half_open], [0, 0], 'k-', linewidth=6)
    ax.plot([half_open, mask_w/2], [0, 0], 'k-', linewidth=6)
    ax.set_xlabel("Width (µm)")
    ax.set_ylabel("Depth (µm)")
    ax.set_title(f"{'異向性 (54.7°)' if mode.startswith('異向性') else '等向性'} 蝕刻剖面示意")
    ax.axhline(0, color='gray', linewidth=0.5)
    ax.invert_yaxis()
    ax.set_aspect('equal', adjustable='datalim')
    st.pyplot(fig)

# =========================================================
# ⑥ Mask 材料建議
# =========================================================
elif page == "⑥ Mask 材料建議":
    st.header("Mask 材料建議工具")
    st.caption("依據所選蝕刻劑，從資料庫中找出「幾乎不被蝕刻」的材料作為建議 Mask")

    etch_choice = st.selectbox("選擇蝕刻劑/條件",
                                etchant_material_rate["Etchant"].tolist())
    row = etchant_material_rate[
        etchant_material_rate["Etchant"] == etch_choice].iloc[0]

    st.write(f"**蝕刻類型：** {row['Type']}　**原設計蝕刻目標：** {row['TargetMaterial']}")

    rates = {m: row[m] for m in MATERIAL_COLS if not pd.isna(row[m])}
    if len(rates) == 0:
        st.warning("此蝕刻劑無材料資料可供比較。")
    else:
        sorted_rates = sorted(rates.items(), key=lambda x: x[1])
        st.subheader("各材料抗蝕刻能力排序（速率越低越適合當 Mask）")
        rank_df = pd.DataFrame(sorted_rates,
                                columns=["材料", "蝕刻速率 (nm/min)"])
        st.dataframe(rank_df, use_container_width=True)

        best_mask = sorted_rates[0]
        st.success(f"✅ 建議 Mask 材料：**{best_mask[0]}** "
                    f"（蝕刻速率僅 {best_mask[1]} nm/min）")

    st.markdown("---")
    st.subheader("硼掺杂 Etch-Stop 定量估算")
    st.latex(r"\frac{R(N_B)}{R_0} \propto \left(\frac{N_{B,ref}}{N_B}\right)^4 \quad (N_B > 2\times10^{19}\,cm^{-3})")
    st.caption("出處：Seidel, H. et al. (1990). \"Anisotropic etching of "
               "crystalline silicon in alkaline solution. Part II: "
               "Influence of dopants.\" J. Electrochem. Soc. 137, 3626-3632。"
               "硼濃度超過約 2×10¹⁹ cm⁻³ 時，蝕刻速率的下降與硼濃度的四次方成反比。")
    NB_REF = 2e19
    NB = st.number_input("硼掺杂濃度 N_B (cm⁻³)", min_value=1e17, max_value=1e22,
                          value=5e19, step=1e19, format="%.3e")
    if NB > NB_REF:
        reduction = (NB_REF / NB) ** 4
        st.info(f"相對於未掺杂矽，蝕刻速率預估降至原速率的 "
                f"**{reduction*100:.4f}%**（僅為定性趨勢估算，非精確定量模型）")
    else:
        st.info("硼濃度低於臨界值 (2×10¹⁹ cm⁻³)，尚不足以產生明顯的 etch-stop 效應。")

# =========================================================
# ⑦ 參考文獻與資料來源
# =========================================================
elif page == "⑦ 參考文獻與資料來源":
    st.header("參考文獻與資料來源")
    st.markdown("""
本應用程式所有動力學參數與公式均可追溯至以下文獻，方便查證：

### 課堂資料
- 課堂 PDF《Etching》。原始表格引用：
  Kovacs, G.T.A., Maluf, N.I., Peterson, K.E. (1998). "Bulk micromachining
  of silicon." *Proceedings of the IEEE*, 86(8), 1536-1551.
  Williams, K. & Muller, R. (作者與年份見 PDF 原文)。

### KOH 動力學參數
- Seidel, H., Csepregi, L., Heuberger, A., Baumgärtel, H. (1990).
  "Anisotropic etching of crystalline silicon in alkaline solution.
  Part I: Orientation dependence and behavior of passivation layer."
  *J. Electrochem. Soc.*, 137, 3612-3626. DOI: 10.1149/1.2086277
- Seidel, H. et al. (1990). "Part II: Influence of dopants."
  *J. Electrochem. Soc.*, 137, 3626-3632.
- Shikida, M., Sato, K., Tokoro, K., Uchikawa, D. (2000). "Differences
  in anisotropic etching properties of KOH and TMAH solutions."
  *Sens. Actuators A*, 80, 179-188.
- Tan, et al. / Wind, et al. -- 引用自二次文獻 (見下方 Handbook)，
  原始期刊出處尚待進一步查證。

### TMAH 動力學參數
- Tabata, O. et al. -- 引用自二次文獻 (見下方 Handbook)，
  原始期刊出處尚待進一步查證。
- Shikida et al. 2000（同上）

### KOH / TMAH / EDP 活化能比較
- Dutta, S., Imran, M., Kumar, P., Pal, R., Datta, P., Chatterjee, R.
  (2011). "Comparison of etch characteristics of KOH, TMAH and EDP
  for bulk micromachining of silicon (110)." *Microsystem Technologies*.

### 二次文獻彙整（Table 22.3, Eq. 22.3, 22.5）
- Gosálvez, M.A., Zubel, I., Viinikka, E. (2015). "Wet Etching of
  Silicon." Chapter 22 in *Handbook of Silicon Based MEMS Materials
  and Technologies* (2nd ed.), Elsevier.
  網址：https://www.sciencedirect.com/science/article/pii/B9780323299657000221

### 乾蝕刻動力學模型
- Ohring, M. (2002). *Materials Science of Thin Films* (2nd ed.),
  Chapter 5 "Plasma and Ion Beam Processing of Thin Films",
  Eq. (5-9a), (5-9b)。

---

### ⚠️ 尚未查證完成的項目（誠實列出）
1. **KOH/TMAH wt%→mol/L 密度換算公式**（Handbook Eq. 22.6 中的
   $a_0=0.9927, a_1=0.8666, a_2=0.3051$ 三個係數）：原文標註見
   Section 22.14，但目前未能取得該章節完整內容，故本版本**未實作
   濃度可調功能**，僅提供固定濃度下的溫度依賴模型。
2. **Tan et al. / Wind et al. / Tabata et al.** 的原始期刊出處：
   目前僅能確認其數據被 Gosálvez et al. (2015) 的 Table 22.3 引用，
   原始期刊、卷期、頁碼尚未逐一查證。
3. **EDP 的前置因子 R0**：文獻 (Dutta et al. 2011) 僅報告活化能
   0.34 eV，未提供 R0 數值，本程式的 R0 是用課堂 PDF 單點資料反推，
   已在程式介面中明確標示。
""")
