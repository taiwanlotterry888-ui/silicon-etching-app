# -*- coding: utf-8 -*-
"""
Silicon Etching Process Evaluator (v3 - 分段聚焦版)
--------------------------------------------------------
本版本依需求「分段開發」，僅聚焦兩個核心功能：
  ② 濕蝕刻動力學計算器（KOH / TMAH / EDP，Arrhenius 模型）
  ③ 乾蝕刻化學動力學計算器（F原子模型）
其餘模組（材料選擇比查詢、異向性剖面模擬、Mask建議）暫時移除，
待這兩個核心功能穩定後再加回來。

本版本修正的問題：
 1. matplotlib 圖表文字全部改為英文 -- Streamlit Cloud 伺服器的
    matplotlib 預設字體 (DejaVu Sans) 不支援中文，中文字會變成方框
    亂碼，此為環境限制而非程式錯誤，故圖表內文字統一使用英文，
    中文說明留在圖表外的文字區塊。
 2. 資料表欄位名稱：程式內部仍用英文變數名（Ea_eV, R0_um_h...）
    方便程式邏輯撰寫，但顯示給使用者的表格會轉成中文/易讀格式，
    不會直接把底線命名的變數名稱顯示出來。
 3. EDP 現在也會出現在「完整資料表」中，並清楚標註它的 R0 是反推值。
 4. 新增說明：HF:HNO3:CH3COOH（HNA系統）目前未納入動力學計算器，
    因為其蝕刻機制為酸性氧化-溶解反應，速率同時取決於 HF/HNO3 的
    比例（而非單一濃度軸），且本次文獻搜索未查到公開發表的 Ea/R0
    數值，故先誠實排除，待查到文獻後再補上。

資料來源總覽（詳見頁面④「參考文獻」）：
 - Seidel, H. et al. (1990). J. Electrochem. Soc. 137, 3612 / 3626
 - Shikida, M. et al. (2000). Sens. Actuators A 80, 179-188
 - Tabata, O. et al.; Tan et al.; Wind et al.
   （以上經 Gosálvez et al. 2015 Table 22.3 彙整引用）
 - Dutta, S. et al. (2011). Microsystem Technologies (KOH/TMAH/EDP活化能比較)
 - Gosálvez, M.A., Zubel, I., Viinikka, E. (2015). "Wet Etching of Silicon."
   Handbook of Silicon Based MEMS Materials and Technologies (2nd ed.), Ch.22
 - Ohring, M. (2002). Materials Science of Thin Films (2nd ed.), Ch.5,
   Eq.(5-9a)(5-9b) -- 乾蝕刻氟原子化學反應模型
 - 課堂 PDF《Etching》(Kovacs, Maluf & Peterson 1998; Williams & Muller)
   -- EDP 單點溫度/速率資料，用於反推 R0

部署方式：
1. 把整個 app.py 覆蓋上傳到 GitHub repo
2. Streamlit Community Cloud 會自動偵測更新並重新部署
"""

import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

st.set_page_config(page_title="Silicon Etching Kinetics Calculator",
                    layout="wide")

KB_EV = 8.617333262e-5  # Boltzmann constant, eV/K

# =========================================================
# 1. 濕蝕刻動力學資料庫 (Arrhenius: R = R0 * exp(-Ea/kB T))
#    R0 以 µm/h 儲存（多數文獻的原始單位），計算時換算成 µm/min
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
EDP_EA = 0.34  # eV, Dutta et al. 2011 (Si(110))
EDP_REF_T_C = 115.0         # PDF General Characteristics 表：EDP 溫度
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


# 建立「完整資料表」用的合併表（含 EDP 反推列），僅供顯示
def build_wet_kinetics_display():
    disp = wet_kinetics.copy()
    edp_row = pd.DataFrame([{
        "Etchant": "EDP", "Concentration": "（未指定，PDF未載明）",
        "Plane": "（未指定）", "Ea_eV": EDP_EA,
        "R0_um_h": edp_r0_um_h(),
        "Source": "Ea: Dutta et al. 2011／R0: 由課堂PDF單點反推（非文獻值）"
    }])
    disp = pd.concat([disp, edp_row], ignore_index=True)
    disp = disp.rename(columns={
        "Etchant": "蝕刻劑",
        "Concentration": "濃度",
        "Plane": "晶面",
        "Ea_eV": "活化能 Ea (eV)",
        "R0_um_h": "前置因子 R0 (µm/h)",
        "Source": "資料來源",
    })
    return disp


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
# Sidebar 導覽
# =========================================================
st.sidebar.title("🔬 Silicon Etching Kinetics")
page = st.sidebar.radio(
    "選擇功能模組",
    [
        "① 首頁說明",
        "② 濕蝕刻動力學計算器（KOH / TMAH / EDP）",
        "③ 乾蝕刻化學動力學計算器（F原子模型）",
        "④ 參考文獻與資料來源",
    ],
)

# =========================================================
# ① 首頁
# =========================================================
if page == "① 首頁說明":
    st.title("矽蝕刻動力學計算器（分段聚焦版）")
    st.markdown("""
本版本聚焦於兩個核心功能，先把物理模型做扎實：

| 模組 | 物理模型 | 可調參數 |
|---|---|---|
| 濕蝕刻動力學計算器 | Arrhenius $R=R_0e^{-E_a/k_BT}$（固定濃度） | 蝕刻劑、晶面、文獻來源、溫度 |
| 乾蝕刻化學動力學計算器 | Arrhenius + 氟原子濃度模型 | 材料、溫度、氟原子濃度 |

### 已知限制（誠實聲明）
- KOH / TMAH 的「濃度」目前不是可連續調整的變數，僅能在文獻報告過的
  固定濃度中選擇（例如 KOH 34/35/37 wt%），因為 wt%→mol/L 換算所需的
  密度公式係數目前查無可靠原始出處。
- EDP 的前置因子 $R_0$ 是由課堂 PDF 單點資料反推，非文獻直接數值。
- **HF:HNO3:CH3COOH（HNA系統）尚未納入本計算器**：這是酸性氧化-溶解
  機制，蝕刻速率同時取決於 HF 與 HNO3 的比例（而非單一濃度軸），
  本次文獻搜索未查到公開發表的 Arrhenius 參數，待查到後再補上。

### 之後會加回來的功能（暫時移除，之後分段補上）
- 蝕刻劑材料選擇比查詢（含 HNA 系統的靜態查表）
- 異向性蝕刻剖面模擬器
- Mask 材料建議工具
""")

# =========================================================
# ② 濕蝕刻動力學計算器
# =========================================================
elif page == "② 濕蝕刻動力學計算器（KOH / TMAH / EDP）":
    st.header("濕蝕刻動力學計算器")
    st.latex(r"R = R_0 \cdot e^{-E_a / k_B T}")
    st.caption("固定濃度模型：濃度已內建於所選文獻資料組中，僅溫度可調整。")

    with st.expander("⚠️ 為什麼沒有 HF:HNO3:CH3COOH（HNA系統）？"):
        st.markdown("""
HNA（HF:HNO3:CH3COOH）是等向性酸性蝕刻，機制是硝酸先氧化矽表面、
氫氟酸再溶解氧化層，速率同時取決於 **HF 與 HNO3 的比例**（业界稱為
"etching triangle"），不是單一濃度軸能描述的簡單系統，也不是像
KOH/TMAH/EDP 這種可以套用單一 Arrhenius 公式的反應。目前文獻搜索
沒有找到這個系統公開發表的 Ea/R0 數值，所以先誠實排除，暫不納入
本計算器，待之後查到文獻參數再補上。
""")

    etchant_choice = st.selectbox("選擇蝕刻劑", ["KOH", "TMAH", "EDP"])

    if etchant_choice in ["KOH", "TMAH"]:
        sub = wet_kinetics[wet_kinetics["Etchant"] == etchant_choice]
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
        st.caption(f"資料來源：{row['Source']}（詳見頁面④）")

        T_C = st.slider("蝕刻溫度 (°C)", 20.0, 130.0, 80.0, 1.0)
        rate = arrhenius_rate_um_min(row["Ea_eV"], row["R0_um_h"], T_C)
        st.success(f"預估蝕刻速率：**{rate:.3f} µm/min**"
                    f"（{rate*60:.2f} µm/h）　＠ {T_C:.0f}°C")

        # Arrhenius 曲線圖（圖表內文字統一用英文，避免中文字體亂碼）
        T_range = np.linspace(20, 130, 200)
        R_range = [arrhenius_rate_um_min(row["Ea_eV"], row["R0_um_h"], t)
                   for t in T_range]
        fig, ax = plt.subplots(figsize=(6, 3.5))
        ax.plot(T_range, R_range, 'b-')
        ax.axvline(T_C, color='red', linestyle='--', alpha=0.6)
        ax.plot([T_C], [rate], 'ro')
        ax.set_xlabel("Temperature (deg C)")
        ax.set_ylabel("Etch rate (um/min)")
        ax.set_yscale('log')
        ax.set_title(f"{etchant_choice} Arrhenius Plot - {row['Plane']} @ {row['Concentration']}")
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
                label='PDF reference point')
        ax.set_xlabel("Temperature (deg C)")
        ax.set_ylabel("Etch rate (um/min)")
        ax.legend()
        ax.set_title("EDP Arrhenius Plot (R0 back-calculated)")
        ax.grid(alpha=0.3)
        st.pyplot(fig)

        st.markdown("---")
        st.subheader("目標深度 → 預估時間")
        target_depth_um2 = st.number_input("目標蝕刻深度 (µm)", min_value=0.0,
                                            value=100.0, step=10.0, key="edp_depth")
        if rate > 0:
            st.info(f"預估所需時間：約 **{target_depth_um2/rate:.1f} 分鐘**")

    st.markdown("---")
    st.subheader("完整濕蝕刻動力學資料表（KOH / TMAH / EDP）")
    st.dataframe(build_wet_kinetics_display(), use_container_width=True)

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
    ax.set_xlabel("Substrate Temperature (deg C)")
    ax.set_ylabel("Si : SiO2 Selectivity")
    ax.set_title(f"Selectivity vs. Temperature (C_F = 10^{C_F_exp:.1f} atoms/cm3)")
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

    st.markdown("---")
    st.subheader("完整乾蝕刻動力學參數表")
    dry_disp = pd.DataFrame([
        {"材料": "Si", "前置因子 A": DRY_PARAMS["Si"]["A"],
         "活化能 Ea (eV)": DRY_PARAMS["Si"]["Ea_eV"]},
        {"材料": "SiO2", "前置因子 A": DRY_PARAMS["SiO2"]["A"],
         "活化能 Ea (eV)": DRY_PARAMS["SiO2"]["Ea_eV"]},
    ])
    st.dataframe(dry_disp, use_container_width=True)

# =========================================================
# ④ 參考文獻與資料來源
# =========================================================
elif page == "④ 參考文獻與資料來源":
    st.header("參考文獻與資料來源")
    st.markdown("""
本應用程式所有動力學參數與公式均可追溯至以下文獻，方便查證：

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
- Tan et al. / Wind et al. -- 引用自二次文獻（見下方 Handbook），
  原始期刊出處尚待進一步查證。

### TMAH 動力學參數
- Tabata, O. et al. -- 引用自二次文獻（見下方 Handbook），
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

### 課堂資料
- 課堂 PDF《Etching》。原始表格引用：
  Kovacs, G.T.A., Maluf, N.I., Peterson, K.E. (1998). "Bulk micromachining
  of silicon." *Proceedings of the IEEE*, 86(8), 1536-1551.
  Williams, K. & Muller, R.（作者與年份見 PDF 原文）。
  -- EDP 單點溫度/速率資料 (115°C, 0.75 µm/min)，用於反推 R0。

---

### ⚠️ 尚未查證完成的項目（誠實列出）
1. **KOH/TMAH wt%→mol/L 密度換算公式**（$a_0=0.9927, a_1=0.8666,
   a_2=0.3051$ 三個係數）：原文標註見 Handbook Section 22.14，
   但目前未能取得該章節完整內容，故本版本未實作濃度可調功能。
2. **Tan et al. / Wind et al. / Tabata et al.** 的原始期刊出處：
   目前僅能確認其數據被 Gosálvez et al. (2015) 的 Table 22.3 引用，
   原始期刊、卷期、頁碼尚未逐一查證。
3. **EDP 的前置因子 R0**：文獻 (Dutta et al. 2011) 僅報告活化能
   0.34 eV，未提供 R0 數值，本程式的 R0 是用課堂 PDF 單點資料反推。
4. **HF:HNO3:CH3COOH（HNA系統）**：尚未查到公開發表的 Arrhenius
   參數 (Ea, R0)，目前完全未納入本計算器。
""")
