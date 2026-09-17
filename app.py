# -*- coding: utf-8 -*-
"""
Silicon Etching Process Evaluator (v4 - 完整功能版)
--------------------------------------------------------
在 v3（分段聚焦版：濕蝕刻動力學計算器 + 乾蝕刻化學動力學計算器）
的基礎上，補回原始藍圖中規劃、之後要加回來的四個模組：

  ④ 蝕刻劑 × 材料選擇比查詢（雙向查詢）
  ⑤ 異向性蝕刻剖面模擬器（54.7° 幾何計算，matplotlib 繪圖）
  ③ 乾蝕刻/電漿計算器（v3 已有，本版保留不動）
  ⑥ Mask 材料建議 + 蝕刻時間反算

v3 既有的三項修正原封不動保留：
 1. matplotlib 圖表文字全部使用英文（Streamlit Cloud 預設字體不支援中文）。
 2. 顯示用資料表欄位為中文，程式內部變數維持英文命名。
 3. EDP 加入完整資料表，並清楚標註 R0 為 PDF 單點反推值。
 4. HF:HNO3:CH3COOH（HNA）未納入 Arrhenius 動力學計算器的說明保留，
    但本版把它的「操作點」數據（非活化能模型）納入新增的④和⑥模組，
    因為 PDF 的 General Characteristics 表格本來就是用單一操作點
    （溫度、蝕刻速率）描述 HNA，而不是 Arrhenius 參數，這兩者不衝突。

新增資料來源：
 - PDF《Etching》p.2 表格 "General Characteristics of Silicon Etching
   Operations"（原始出處：Kovacs, Maluf & Peterson (1998), Proceedings
   of the IEEE, 86(8), 1536-1551）
   → 用於：④材料選擇比查詢（製程層級）、⑥Mask建議與蝕刻時間反算
 - PDF《Etching》p.7 表格 "Comparison of Etch Rates for Selected
   Etchants and Target Materials"（原始出處：Williams, K. & Muller, R.）
   → 用於：④蝕刻劑 × 材料選擇比查詢（材料層級）
 - 54.7° 幾何：PDF《Etching》"Anisotropic Etching" 頁的示意圖
   （{111} 面與 {100} 面夾角 54.7°，KOH/TMAH 沿 {111} 面停止的異向性
   蝕刻幾何）→ 用於：⑤異向性蝕刻剖面模擬器

⚠️ 誠實聲明：④和⑥模組使用的 PDF 表格數值是「單一操作點」的參考值
（例如某個溫度下量測到的速率），不是 Arrhenius 的 (Ea, R0) 兩參數，
所以這兩個模組**不能**像②的計算器一樣連續調整溫度來外插速率，只能
在 PDF 報告的操作條件附近使用，程式內以文字提醒使用者這個限制。

部署方式：
1. 把整個 app.py 覆蓋上傳到 GitHub repo
2. Streamlit Community Cloud 會自動偵測更新並重新部署
"""

import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, NullLocator

st.set_page_config(page_title="Silicon Etching Kinetics Calculator",
                    layout="wide")

KB_EV = 8.617333262e-5  # Boltzmann constant, eV/K
COT_547 = 1.0 / np.tan(np.deg2rad(54.7))  # {111}/{100} 夾角的餘切值

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
EDP_EXTRAPOLATION_WARN_C = 15.0  # 距參考點超過此溫差 (°C) 時另外提出外插警告


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
# 3. 製程層級「操作點」資料表（PDF p.2 General Characteristics）
#    這些是單一操作條件下的參考值，不是 Arrhenius 參數。
#    數值取表格區間的代表值（中位數/上界），供選擇比估算與
#    Mask 蝕刻時間反算使用，皆在 UI 上標註「近似值」。
# =========================================================
general_char = pd.DataFrame([
    {"Process": "HF:HNO3:CH3COOH (HNA)", "Type": "Wet",
     "T_C_range": "25", "EtchRate_um_min_range": "1-20",
     "Selectivity_111_100": "—",
     "Nitride_rate_nm_min_range": "Low", "Nitride_rate_nm_min_point": 1.0,
     "SiO2_rate_nm_min_range": "10-30", "SiO2_rate_nm_min_point": 20.0,
     "pplus_etch_stop": "No"},
    {"Process": "KOH", "Type": "Wet",
     "T_C_range": "70-90", "EtchRate_um_min_range": "0.5-2",
     "Selectivity_111_100": "100:1",
     "Nitride_rate_nm_min_range": "<1", "Nitride_rate_nm_min_point": 0.5,
     "SiO2_rate_nm_min_range": "10", "SiO2_rate_nm_min_point": 10.0,
     "pplus_etch_stop": "Yes"},
    {"Process": "EDP (ethylene-diamine pyrocatechol)", "Type": "Wet",
     "T_C_range": "115", "EtchRate_um_min_range": "0.75",
     "Selectivity_111_100": "35:1",
     "Nitride_rate_nm_min_range": "0.1", "Nitride_rate_nm_min_point": 0.1,
     "SiO2_rate_nm_min_range": "0.2", "SiO2_rate_nm_min_point": 0.2,
     "pplus_etch_stop": "Yes"},
    {"Process": "N(CH3)4OH (TMAH)", "Type": "Wet",
     "T_C_range": "90", "EtchRate_um_min_range": "0.5-1.5",
     "Selectivity_111_100": "50:1",
     "Nitride_rate_nm_min_range": "<0.1", "Nitride_rate_nm_min_point": 0.05,
     "SiO2_rate_nm_min_range": "<0.1", "SiO2_rate_nm_min_point": 0.05,
     "pplus_etch_stop": "Yes"},
    {"Process": "SF6 (plasma)", "Type": "Dry",
     "T_C_range": "0-100", "EtchRate_um_min_range": "0.1-0.5",
     "Selectivity_111_100": "—",
     "Nitride_rate_nm_min_range": "200", "Nitride_rate_nm_min_point": 200.0,
     "SiO2_rate_nm_min_range": "10", "SiO2_rate_nm_min_point": 10.0,
     "pplus_etch_stop": "No"},
    {"Process": "SF6/C4F8 (DRIE)", "Type": "Dry",
     "T_C_range": "20-80", "EtchRate_um_min_range": "1-3",
     "Selectivity_111_100": "—",
     "Nitride_rate_nm_min_range": "200", "Nitride_rate_nm_min_point": 200.0,
     "SiO2_rate_nm_min_range": "10", "SiO2_rate_nm_min_point": 10.0,
     "pplus_etch_stop": "No"},
])

GENERAL_CHAR_SOURCE = ("Kovacs, G.T.A., Maluf, N.I., Peterson, K.E. (1998). "
                        "\"Bulk micromachining of silicon.\" Proceedings of "
                        "the IEEE, 86(8), 1536-1551（課堂 PDF 轉引）")


def build_general_char_display():
    disp = general_char[[
        "Process", "Type", "T_C_range", "EtchRate_um_min_range",
        "Selectivity_111_100", "Nitride_rate_nm_min_range",
        "SiO2_rate_nm_min_range", "pplus_etch_stop"
    ]].rename(columns={
        "Process": "製程", "Type": "濕/乾",
        "T_C_range": "溫度 (°C)", "EtchRate_um_min_range": "Si 蝕刻速率 (µm/min)",
        "Selectivity_111_100": "{111}/{100} 選擇比",
        "Nitride_rate_nm_min_range": "氮化矽蝕刻速率 (nm/min)",
        "SiO2_rate_nm_min_range": "SiO2 蝕刻速率 (nm/min)",
        "pplus_etch_stop": "p++ 蝕刻終止",
    })
    return disp


# =========================================================
# 4. 蝕刻劑 × 目標材料選擇比資料表（PDF p.7 Isotropic Etchants 頁）
#    來源：Williams, K. & Muller, R.（課堂 PDF 轉引），單位 nm/min。
#    ">X" 取為下界估計值 X（實際更快）；"<X" 取為上界估計值 X（實際更慢）；
#    "—" 表示原表未報告，以 NaN 處理。這些皆為「單次量測操作點」，
#    非隨溫度連續變化的模型，UI 會標註此限制。
# =========================================================
material_etch = pd.DataFrame([
    {"Etchant": "Concentrated HF (49%)", "Category": "Wet", "DesignedFor": "Silicon oxides",
     "Polysilicon_n+": 0, "Polysilicon_undoped": 0, "SiO2": 2300, "SiNitride": 14,
     "PSG_annealed": 3600, "Aluminum": 4.2, "Titanium": 1000, "Photoresist": 0,
     "Note": ""},
    {"Etchant": "25:1 HF:H2O", "Category": "Wet", "DesignedFor": "Silicon oxides",
     "Polysilicon_n+": 0, "Polysilicon_undoped": 0, "SiO2": 9.7, "SiNitride": 0.6,
     "PSG_annealed": 150, "Aluminum": np.nan, "Titanium": np.nan, "Photoresist": 0,
     "Note": ""},
    {"Etchant": "5:1 BHF", "Category": "Wet", "DesignedFor": "Silicon oxides",
     "Polysilicon_n+": 9, "Polysilicon_undoped": 2, "SiO2": 100, "SiNitride": 0.9,
     "PSG_annealed": 440, "Aluminum": 140, "Titanium": 1000, "Photoresist": 0,
     "Note": "Titanium 為 >1000 下界估計"},
    {"Etchant": "Silicon etchant (126HNO3:60H2O:5NH4F)", "Category": "Wet", "DesignedFor": "Silicon",
     "Polysilicon_n+": 310, "Polysilicon_undoped": 100, "SiO2": 9, "SiNitride": 0.2,
     "PSG_annealed": 170, "Aluminum": 400, "Titanium": 300, "Photoresist": 0,
     "Note": ""},
    {"Etchant": "Aluminum etchant (16H3PO4:1HNO3:2H2O)", "Category": "Wet", "DesignedFor": "Aluminum",
     "Polysilicon_n+": 1, "Polysilicon_undoped": np.nan, "SiO2": 1, "SiNitride": 0,
     "PSG_annealed": 1, "Aluminum": 660, "Titanium": 0, "Photoresist": 0,
     "Note": "Polysilicon/SiO2/PSG 為 <1 上界估計"},
    {"Etchant": "Titanium etchant (20H2O:1H2O2:1HF)", "Category": "Wet", "DesignedFor": "Titanium",
     "Polysilicon_n+": 1.2, "Polysilicon_undoped": np.nan, "SiO2": 12, "SiNitride": 0,
     "PSG_annealed": 210, "Aluminum": 10, "Titanium": 880, "Photoresist": 0,
     "Note": "Aluminum 為 >10 下界估計"},
    {"Etchant": "Piranha (50H2SO4:1H2O2)", "Category": "Wet", "DesignedFor": "Metal/organic cleaning",
     "Polysilicon_n+": 0, "Polysilicon_undoped": 0, "SiO2": 0, "SiNitride": 0,
     "PSG_annealed": 0, "Aluminum": 10, "Titanium": 240, "Photoresist": 10,
     "Note": "Aluminum/Photoresist 為 >10 下界估計"},
    {"Etchant": "Acetone (CH3COOH)", "Category": "Wet", "DesignedFor": "Photoresist",
     "Polysilicon_n+": 0, "Polysilicon_undoped": 0, "SiO2": 0, "SiNitride": 0,
     "PSG_annealed": 0, "Aluminum": 0, "Titanium": 0, "Photoresist": 4000,
     "Note": "Photoresist 為 >4000 下界估計"},
    {"Etchant": "CF4+CHF3+He, 450W (RIE)", "Category": "Dry", "DesignedFor": "Silicon oxides",
     "Polysilicon_n+": 190, "Polysilicon_undoped": 210, "SiO2": 470, "SiNitride": 180,
     "PSG_annealed": 620, "Aluminum": np.nan, "Titanium": 1000, "Photoresist": 220,
     "Note": ""},
    {"Etchant": "SF6+He, 100W (RIE)", "Category": "Dry", "DesignedFor": "Silicon nitrides",
     "Polysilicon_n+": 73, "Polysilicon_undoped": 67, "SiO2": 31, "SiNitride": 32,
     "PSG_annealed": 61, "Aluminum": np.nan, "Titanium": 1000, "Photoresist": 69,
     "Note": ""},
    {"Etchant": "SF6, 12.5W (RIE)", "Category": "Dry", "DesignedFor": "Thin silicon nitrides",
     "Polysilicon_n+": 170, "Polysilicon_undoped": 280, "SiO2": 110, "SiNitride": 280,
     "PSG_annealed": 140, "Aluminum": np.nan, "Titanium": 1000, "Photoresist": 310,
     "Note": ""},
    {"Etchant": "O2, 400W (plasma ashing)", "Category": "Dry", "DesignedFor": "Photoresist ashing",
     "Polysilicon_n+": 0, "Polysilicon_undoped": 0, "SiO2": 0, "SiNitride": 0,
     "PSG_annealed": 0, "Aluminum": 0, "Titanium": 0, "Photoresist": 340,
     "Note": ""},
])

MATERIAL_ETCH_SOURCE = "Williams, K. & Muller, R.（課堂 PDF《Etching》Isotropic Etchants 頁轉引）"
MATERIAL_COLS = ["Polysilicon_n+", "Polysilicon_undoped", "SiO2", "SiNitride",
                  "PSG_annealed", "Aluminum", "Titanium", "Photoresist"]
MATERIAL_COLS_ZH = {
    "Polysilicon_n+": "多晶矽 (n+摻雜)", "Polysilicon_undoped": "多晶矽 (未摻雜)",
    "SiO2": "二氧化矽", "SiNitride": "氮化矽", "PSG_annealed": "磷矽玻璃(退火後)",
    "Aluminum": "鋁", "Titanium": "鈦", "Photoresist": "光阻(OCG-820PR)",
}


def build_material_etch_display():
    disp = material_etch.rename(columns={
        "Etchant": "蝕刻劑", "Category": "濕/乾", "DesignedFor": "主要設計目標材料",
        "Note": "備註（估計值說明）", **MATERIAL_COLS_ZH,
    })
    return disp


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
        "④ 蝕刻劑 × 材料選擇比查詢",
        "⑤ 異向性蝕刻剖面模擬器（54.7°幾何）",
        "⑥ Mask材料建議與蝕刻時間反算",
        "⑦ 參考文獻與資料來源",
    ],
)

# =========================================================
# ① 首頁
# =========================================================
if page == "① 首頁說明":
    st.title("矽蝕刻動力學計算器（完整功能版）")
    st.markdown("""
本應用程式協助評估矽蝕刻製程的蝕刻速率、材料選擇比與蝕刻剖面，
分為七個模組：

| 模組 | 物理模型 / 資料來源 | 可調參數 |
|---|---|---|
| ② 濕蝕刻動力學計算器 | Arrhenius $R=R_0e^{-E_a/k_BT}$（KOH/TMAH/EDP，固定濃度） | 蝕刻劑、晶面、文獻來源、溫度 |
| ③ 乾蝕刻化學動力學計算器 | Arrhenius + 氟原子濃度模型（Ohring 2002） | 材料、溫度、氟原子濃度 |
| ④ 蝕刻劑 × 材料選擇比查詢 | PDF 操作點資料表（Williams & Muller） | 雙向查詢：選蝕刻劑看材料、或選材料看蝕刻劑 |
| ⑤ 異向性蝕刻剖面模擬器 | 54.7° {111}/{100} 幾何 | 遮罩開口寬度、蝕刻深度（或時間） |
| ⑥ Mask材料建議與蝕刻時間反算 | PDF 操作點資料表 + ②的 Arrhenius 計算器 | 目標深度 ↔ 所需時間 ↔ Mask 厚度 |

### 已知限制（誠實聲明）
- KOH / TMAH 的「濃度」目前不是可連續調整的變數，僅能在文獻報告過的
  固定濃度中選擇（例如 KOH 34/35/37 wt%），因為 wt%→mol/L 換算所需的
  密度公式係數目前查無可靠原始出處。
- EDP 的前置因子 $R_0$ 是由課堂 PDF 單點資料反推，非文獻直接數值。
- **HF:HNO3:CH3COOH（HNA系統）沒有 Arrhenius 動力學模型**：這是酸性
  氧化-溶解機制，蝕刻速率同時取決於 HF 與 HNO3 的比例（而非單一濃度
  軸），本次文獻搜索未查到公開發表的 Arrhenius 參數。不過④和⑥模組
  仍會用 PDF 表格中 HNA 的「單一操作點」數值（25°C 下的速率）做選擇比
  估算，因為那本來就是操作點資料，不需要 Arrhenius 模型。
- ④和⑥模組使用的 PDF 表格數值都是**單一操作條件下的參考值**，不能
  像②一樣連續外插溫度；含 "<" / ">" 的數值一律以近似值處理，並在表格
  中加註。
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
氫氟酸再溶解氧化層，速率同時取決於 **HF 與 HNO3 的比例**（業界稱為
"etching triangle"），不是單一濃度軸能描述的簡單系統，也不是像
KOH/TMAH/EDP 這種可以套用單一 Arrhenius 公式的反應。目前文獻搜索
沒有找到這個系統公開發表的 Ea/R0 數值，所以先誠實排除，暫不納入
本計算器，待之後查到文獻參數再補上。它在 PDF 表格中的單一操作點
數值（25°C, 1-20 µm/min）仍可在④和⑥模組中用於選擇比估算。
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
        st.caption(f"資料來源：{row['Source']}（詳見頁面⑦）")

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
        ax.set_xlabel("Temperature (\u00b0C)")
        ax.set_ylabel("Etch rate (\u00b5m/min)")
        ax.set_yscale('log')
        # 只顯示 10 的整數次方刻度（10^-1, 10^0, 10^1...），關閉次刻度，
        # 避免 log 座標軸出現一堆沒有標籤、看起來參差不齊的細刻度。
        ax.yaxis.set_major_locator(LogLocator(base=10.0))
        ax.yaxis.set_minor_locator(NullLocator())
        ax.set_title(f"{etchant_choice} Arrhenius Plot ({row['Plane']} plane, {row['Concentration']})")
        ax.grid(alpha=0.3, which='major')
        st.pyplot(fig)

        # ---- {100}/{111} 選擇比：與講義「100:1」典型值比較 ----
        same_group = wet_kinetics[
            (wet_kinetics["Etchant"] == etchant_choice) &
            (wet_kinetics["Concentration"] == row["Concentration"]) &
            (wet_kinetics["Source"] == row["Source"])
        ]
        planes_available = set(same_group["Plane"])
        if {"{100}", "{111}"}.issubset(planes_available):
            row100 = same_group[same_group["Plane"] == "{100}"].iloc[0]
            row111 = same_group[same_group["Plane"] == "{111}"].iloc[0]
            rate100 = arrhenius_rate_um_min(row100["Ea_eV"], row100["R0_um_h"], T_C)
            rate111 = arrhenius_rate_um_min(row111["Ea_eV"], row111["R0_um_h"], T_C)
            if rate111 > 0:
                selectivity = rate100 / rate111
                st.markdown("---")
                st.subheader("{100} / {111} 選擇比（與講義典型值比較）")
                col1, col2, col3 = st.columns(3)
                with col1:
                    st.metric("{100} 蝕刻速率", f"{rate100:.3f} µm/min")
                with col2:
                    st.metric("{111} 蝕刻速率", f"{rate111:.4f} µm/min")
                with col3:
                    st.metric("計算出的選擇比", f"{selectivity:.1f} : 1")

                with st.expander("💡 為什麼跟講義標示的「100:1」不一樣？"):
                    st.markdown(f"""
講義 General Characteristics 表格中，KOH 的 {{111}}/{{100}} 選擇比標示
為 **100:1**，但這是**橫跨多種操作條件（不同濃度、不同溫度）的
「數量級與典型最大值概估」**，用意是讓人快速掌握「KOH 對 {{111}} 面
的蝕刻明顯比 {{100}} 面慢很多」這個定性趨勢，不是針對某一個特定濃度
的精確數字。

本計算器套用的則是 **{row['Source']}** 在 **{row['Concentration']}、
{T_C:.0f}°C** 這個具體條件下的 Arrhenius 精確經驗公式，兩者本來就是
回答不同層次的問題：

- **濃度是關鍵變數**：KOH 的 {{111}}/{{100}} 選擇比對溶液濃度非常
  敏感。文獻報告中，較低濃度（例如 20 wt% 左右）選擇比可以達到
  100:1 甚至更高；但在較高濃度（例如本組所用的 {row['Concentration']}）
  下，選擇比通常會下降，落在約 50:1～70:1 之間。
- 本次算出的 **{selectivity:.1f} : 1**，正好落在 34 wt% KOH 文獻報告的
  合理實驗範圍內，跟講義的「100:1」**並不矛盾**——講義給的是跨條件
  的概估值，這裡給的是特定濃度/溫度下的精確計算值。
- 換句話說：兩個數字不同不代表任一邊算錯，而是「數量級概估」與
  「特定條件精確模型」性質不同，使用時請注意不要直接互相取代或
  拿來互相「校正」。
""")
        else:
            st.caption("（此文獻資料組沒有同時提供 {100} 與 {111} 的資料，"
                       "無法在此計算選擇比。）")

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

        temp_diff = abs(T_C - EDP_REF_T_C)
        if temp_diff > EDP_EXTRAPOLATION_WARN_C:
            st.warning(f"""
⚠️ **外插警告**：目前溫度（{T_C:.0f}°C）與反推 R0 所用的參考點
（{EDP_REF_T_C:.0f}°C）相差 **{temp_diff:.0f}°C**。

R0 是只用**單一個溫度、單一個速率數值**反推出來的，並不是像 KOH/TMAH
那樣由多筆不同溫度的實驗數據回歸得到 —— 換句話說，Arrhenius 曲線的
「斜率」（活化能 Ea）雖然是文獻值，但曲線的「高度」（R0）只在
{EDP_REF_T_C:.0f}°C 這一點被驗證過。離這個參考點越遠，計算出的速率
就越接近**單純外插**，沒有實驗數據可以佐證其準確度，僅供粗略參考，
不建議直接用於製程設計。
""")
        elif temp_diff > 0:
            st.info(f"目前溫度與反推 R0 所用的參考點（{EDP_REF_T_C:.0f}°C）"
                    f"相差 {temp_diff:.0f}°C，仍在鄰近範圍內，但請記得 R0 "
                    "本身只由單一資料點反推，準確度仍有限。")

        T_range = np.linspace(60, 130, 200)
        R_range = [arrhenius_rate_um_min(EDP_EA, R0, t) for t in T_range]
        fig, ax = plt.subplots(figsize=(6, 3.5))
        ax.plot(T_range, R_range, 'g-')
        ax.axvline(T_C, color='red', linestyle='--', alpha=0.6)
        ax.plot([T_C], [rate], 'ro')
        ax.plot([EDP_REF_T_C], [EDP_REF_RATE_UM_MIN], 'k*', markersize=12,
                label='PDF reference point')
        ax.set_xlabel("Temperature (\u00b0C)")
        ax.set_ylabel("Etch rate (\u00b5m/min)")
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
    ax.set_xlabel("Substrate Temperature (\u00b0C)")
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
# ④ 蝕刻劑 × 材料選擇比查詢（雙向）
# =========================================================
elif page == "④ 蝕刻劑 × 材料選擇比查詢":
    st.header("蝕刻劑 × 材料選擇比查詢")
    st.caption("資料來源：" + MATERIAL_ETCH_SOURCE + "。單位 nm/min，"
               "皆為單一操作條件下的參考值（非 Arrhenius 模型），"
               "含 '<'/'>' 者以估計值處理，詳見各列「備註」欄。")

    mode = st.radio("查詢方向", ["已知蝕刻劑 → 查各材料蝕刻速率",
                                  "已知目標材料 → 查各蝕刻劑蝕刻速率"],
                     horizontal=True)

    if mode == "已知蝕刻劑 → 查各材料蝕刻速率":
        etchant_pick = st.selectbox("選擇蝕刻劑", material_etch["Etchant"].tolist())
        row = material_etch[material_etch["Etchant"] == etchant_pick].iloc[0]
        st.caption(f"主要設計目標材料：**{row['DesignedFor']}**"
                   + (f"　｜　備註：{row['Note']}" if row["Note"] else ""))

        # 圖表 y 軸標籤統一用英文欄位名稱，避免 Streamlit Cloud 中文字體亂碼；
        # 中文材料名稱已在上方 caption 與下方完整資料表中呈現。
        valid_cols = [c for c in MATERIAL_COLS if not pd.isna(row[c])]
        valid_rates = [float(row[c]) for c in valid_cols]

        fig, ax = plt.subplots(figsize=(7, 4))
        ax.barh(valid_cols, valid_rates)
        ax.set_xlabel("Etch rate (nm/min)")
        ax.set_title(f"Etch rate by material - {etchant_pick}")
        ax.set_xscale('symlog')
        ax.grid(alpha=0.3, axis='x')
        st.pyplot(fig)

        st.subheader("與其他材料的選擇比")
        target_mat = st.selectbox("欲蝕刻的目標材料", MATERIAL_COLS,
                                   format_func=lambda c: MATERIAL_COLS_ZH[c])
        protect_mat = st.selectbox("欲保護的材料", MATERIAL_COLS,
                                    format_func=lambda c: MATERIAL_COLS_ZH[c],
                                    index=min(1, len(MATERIAL_COLS)-1))
        r_target = float(row[target_mat]) if not pd.isna(row[target_mat]) else None
        r_protect = float(row[protect_mat]) if not pd.isna(row[protect_mat]) else None
        if r_target is None or r_protect is None:
            st.warning("所選材料在此蝕刻劑的原始表格中為「未報告」，無法計算選擇比。")
        elif r_protect == 0:
            st.success(f"選擇比：**∞ : 1**（{MATERIAL_COLS_ZH[protect_mat]} 幾乎不被此蝕刻劑蝕刻）")
        else:
            st.success(f"{MATERIAL_COLS_ZH[target_mat]} : {MATERIAL_COLS_ZH[protect_mat]} "
                       f"選擇比 ≈ **{r_target/r_protect:.1f} : 1**")

    else:  # 已知材料 → 查蝕刻劑
        material_pick = st.selectbox("選擇目標材料", MATERIAL_COLS,
                                      format_func=lambda c: MATERIAL_COLS_ZH[c])
        sub = material_etch[["Etchant", "Category", material_pick]].copy()
        sub = sub.dropna(subset=[material_pick]).sort_values(
            material_pick, ascending=False)
        sub = sub.rename(columns={"Etchant": "蝕刻劑", "Category": "濕/乾",
                                   material_pick: f"{MATERIAL_COLS_ZH[material_pick]} 蝕刻速率 (nm/min)"})
        st.dataframe(sub, use_container_width=True, hide_index=True)

        if len(sub) > 0:
            fig, ax = plt.subplots(figsize=(7, 4))
            ax.barh(sub["蝕刻劑"], sub.iloc[:, -1])
            ax.set_xlabel("Etch rate (nm/min)")
            ax.set_title(f"Etchant ranking for {material_pick}")
            ax.set_xscale('symlog')
            ax.invert_yaxis()
            ax.grid(alpha=0.3, axis='x')
            st.pyplot(fig)
        st.caption("速率最高者最適合用來「蝕刻」此材料；速率最低（含0）者最適合"
                   "用來當作蝕刻此材料製程的「保護層／Mask」。")

    st.markdown("---")
    st.subheader("完整材料選擇比資料表")
    st.dataframe(build_material_etch_display(), use_container_width=True)

    st.markdown("---")
    st.subheader("製程層級選擇比（PDF General Characteristics 表）")
    st.caption("資料來源：" + GENERAL_CHAR_SOURCE)
    st.dataframe(build_general_char_display(), use_container_width=True)

# =========================================================
# ⑤ 異向性蝕刻剖面模擬器
# =========================================================
elif page == "⑤ 異向性蝕刻剖面模擬器（54.7°幾何）":
    st.header("異向性蝕刻剖面模擬器")
    st.caption("適用於 KOH / TMAH 沿 {100} 晶圓、遮罩對準 <110> 方向的異向性"
               "濕蝕刻。蝕刻側壁停在 {111} 面，與晶圓表面夾角固定為 "
               "**54.7°**（PDF《Etching》Anisotropic Etching 頁示意圖）。")
    st.latex(r"\tan(54.7^\circ) = \frac{\text{depth}}{\text{lateral undercut per side}}")

    depth_mode = st.radio("蝕刻深度輸入方式",
                           ["直接輸入深度", "由②的 Arrhenius 計算器換算時間→深度"],
                           horizontal=True)

    col1, col2 = st.columns(2)
    with col1:
        mask_opening_um = st.number_input("遮罩開口寬度 W (µm)", min_value=1.0,
                                           value=200.0, step=10.0)
    with col2:
        wafer_thickness_um = st.number_input("矽晶圓厚度（供剖面圖比例參考，µm）",
                                              min_value=10.0, value=400.0, step=10.0)

    if depth_mode == "直接輸入深度":
        depth_um = st.number_input("目標蝕刻深度 D (µm)", min_value=0.0,
                                    value=100.0, step=5.0)
        rate_note = None
    else:
        etchant_choice5 = st.selectbox("選擇蝕刻劑（用②的文獻資料組）",
                                        ["KOH", "TMAH"], key="p5_etchant")
        sub5 = wet_kinetics[(wet_kinetics["Etchant"] == etchant_choice5) &
                             (wet_kinetics["Plane"] == "{100}")]
        options5 = sub5.apply(
            lambda r: f"{r['Concentration']} | {r['Source']}", axis=1).tolist()
        choice5 = st.selectbox("選擇 {100} 文獻資料組（決定垂直蝕刻速率）", options5)
        row5 = sub5.iloc[options5.index(choice5)]
        T_C5 = st.slider("蝕刻溫度 (°C)", 20.0, 130.0, 80.0, 1.0, key="p5_temp")
        etch_time_min = st.number_input("蝕刻時間 (分鐘)", min_value=0.0,
                                         value=60.0, step=5.0)
        rate5 = arrhenius_rate_um_min(row5["Ea_eV"], row5["R0_um_h"], T_C5)
        depth_um = rate5 * etch_time_min
        rate_note = (f"{etchant_choice5} @ {row5['Concentration']}, {T_C5:.0f}°C → "
                     f"{{100}} 垂直蝕刻速率 = {rate5:.3f} µm/min")
        st.info(f"{rate_note}　⇒　蝕刻深度 = {depth_um:.1f} µm")

    # ---- 54.7° 幾何計算 ----
    max_depth_um = (mask_opening_um / 2.0) * np.tan(np.deg2rad(54.7))
    is_v_groove = depth_um >= max_depth_um
    if is_v_groove:
        actual_depth_um = max_depth_um
        bottom_width_um = 0.0
    else:
        actual_depth_um = depth_um
        undercut_per_side_um = actual_depth_um * COT_547
        bottom_width_um = mask_opening_um - 2 * undercut_per_side_um

    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("最大可達深度（V型溝槽自限深度）", f"{max_depth_um:.1f} µm")
    with col2:
        st.metric("實際蝕刻深度", f"{actual_depth_um:.1f} µm")
    with col3:
        st.metric("溝槽底部寬度", f"{bottom_width_um:.1f} µm")

    if is_v_groove:
        st.warning(f"輸入的目標深度（{depth_um:.1f} µm）已達到或超過此開口寬度下"
                   f"的自限深度（{max_depth_um:.1f} µm）。異向性蝕刻在兩側 {{111}} "
                   f"面交會後會**自我終止**，實際只能蝕到 V 型溝槽底部，無法再更深。")

    # ---- 剖面圖（英文標籤，避免中文字體亂碼）----
    fig, ax = plt.subplots(figsize=(7, 5))
    half_top = mask_opening_um / 2.0
    if is_v_groove:
        xs = [-half_top - 20, -half_top, 0, half_top, half_top + 20]
        ys = [0, 0, -actual_depth_um, 0, 0]
    else:
        half_bottom = bottom_width_um / 2.0
        xs = [-half_top - 20, -half_top, -half_bottom, half_bottom,
              half_top, half_top + 20]
        ys = [0, 0, -actual_depth_um, -actual_depth_um, 0, 0]
    ax.plot(xs, ys, 'b-', linewidth=2)
    ax.fill_between(xs, ys, min(ys) - 20, color='lightblue', alpha=0.5)

    # 遮罩層
    mask_thickness_draw = max(actual_depth_um * 0.03, 2)
    ax.add_patch(plt.Rectangle((-half_top - 20, 0), 20, mask_thickness_draw,
                                color='gray'))
    ax.add_patch(plt.Rectangle((half_top, 0), 20, mask_thickness_draw,
                                color='gray', label='Mask layer'))

    ax.axhline(0, color='black', linewidth=0.8)
    ax.set_xlabel("Lateral position (um)")
    ax.set_ylabel("Depth (um)")
    ax.set_title(f"Anisotropic Etch Profile (54.7 deg) - W_mask = {mask_opening_um:.0f} um")
    ax.set_ylim(min(-wafer_thickness_um * 0.3, min(ys) - 20), 20)
    ax.set_aspect('equal', adjustable='box')
    ax.grid(alpha=0.3)
    ax.legend(loc='upper right')
    st.pyplot(fig)

    st.markdown("---")
    st.caption("幾何公式：單側側向 undercut = 深度 × cot(54.7°)；"
               "底部寬度 = 遮罩開口寬度 − 2 × 單側 undercut；"
               "自限深度（V型溝槽）= (遮罩開口寬度/2) × tan(54.7°)。"
               "來源：PDF《Etching》Anisotropic Etching 頁示意圖 (b)。")

# =========================================================
# ⑥ Mask材料建議與蝕刻時間反算
# =========================================================
elif page == "⑥ Mask材料建議與蝕刻時間反算":
    st.header("Mask 材料建議與蝕刻時間反算")
    st.caption("結合②的 Arrhenius 蝕刻速率計算器，與 PDF General Characteristics "
               "表格中 SiO2 / 氮化矽的操作點蝕刻速率，估算 Mask 需要多厚才能撐過"
               "整個蝕刻製程。⚠️ Mask 材料的蝕刻速率是單一操作點數值，不會隨溫度"
               "連續調整，計算結果僅供概略設計參考。")

    etchant_choice6 = st.selectbox("選擇蝕刻劑（濕蝕刻）", ["KOH", "TMAH", "EDP"],
                                    key="p6_etchant")

    if etchant_choice6 in ["KOH", "TMAH"]:
        sub6 = wet_kinetics[(wet_kinetics["Etchant"] == etchant_choice6) &
                             (wet_kinetics["Plane"] == "{100}")]
        options6 = sub6.apply(
            lambda r: f"{r['Concentration']} | {r['Source']}", axis=1).tolist()
        choice6 = st.selectbox("選擇 {100} 文獻資料組", options6)
        row6 = sub6.iloc[options6.index(choice6)]
        T_C6 = st.slider("蝕刻溫度 (°C)", 20.0, 130.0, 80.0, 1.0, key="p6_temp")
        si_rate_um_min = arrhenius_rate_um_min(row6["Ea_eV"], row6["R0_um_h"], T_C6)
        process_key = "KOH" if etchant_choice6 == "KOH" else "N(CH3)4OH (TMAH)"
    else:
        st.info(f"EDP 固定使用課堂 PDF 單點條件：{EDP_REF_T_C:.0f}°C。")
        R0_6 = edp_r0_um_h()
        si_rate_um_min = arrhenius_rate_um_min(EDP_EA, R0_6, EDP_REF_T_C)
        process_key = "EDP (ethylene-diamine pyrocatechol)"

    gc_row = general_char[general_char["Process"] == process_key].iloc[0]
    nitride_rate_nm_min = gc_row["Nitride_rate_nm_min_point"]
    sio2_rate_nm_min = gc_row["SiO2_rate_nm_min_point"]

    st.markdown(f"**Si 蝕刻速率：{si_rate_um_min:.3f} µm/min**　｜　"
               f"Mask 蝕刻速率參考點（{process_key} 操作條件）：SiO2 = "
               f"{sio2_rate_nm_min} nm/min，氮化矽 = {nitride_rate_nm_min} nm/min　"
               f"（來源：{GENERAL_CHAR_SOURCE}）")

    safety_factor = st.slider("Mask 厚度安全係數（預留裕度）", 1.0, 3.0, 1.3, 0.1)

    calc_mode = st.radio("反算方向",
                          ["① 給定目標深度 → 求所需時間與Mask厚度",
                           "② 給定現有Mask厚度 → 求最大允許時間與深度"],
                          horizontal=True)

    if calc_mode == "① 給定目標深度 → 求所需時間與Mask厚度":
        target_depth_um6 = st.number_input("目標 Si 蝕刻深度 (µm)", min_value=0.0,
                                            value=100.0, step=10.0, key="p6_depth")
        if si_rate_um_min > 0:
            req_time_min = target_depth_um6 / si_rate_um_min
            sio2_needed_nm = sio2_rate_nm_min * req_time_min * safety_factor
            nitride_needed_nm = nitride_rate_nm_min * req_time_min * safety_factor

            st.success(f"預估所需蝕刻時間：**{req_time_min:.1f} 分鐘**")
            col1, col2 = st.columns(2)
            with col1:
                st.metric("SiO2 Mask 所需厚度（含安全係數）",
                          f"{sio2_needed_nm:.0f} nm")
            with col2:
                st.metric("氮化矽 Mask 所需厚度（含安全係數）",
                          f"{nitride_needed_nm:.0f} nm")

            if nitride_needed_nm < sio2_needed_nm:
                st.info("💡 **建議使用氮化矽 (Si3N4) 當 Mask**："
                       f"在此製程/溫度條件下，氮化矽被蝕刻的速度比 SiO2 慢，"
                       f"所需 Mask 厚度較薄（{nitride_needed_nm:.0f} nm vs. "
                       f"{sio2_needed_nm:.0f} nm）。")
            elif sio2_needed_nm < nitride_needed_nm:
                st.info("💡 **建議使用 SiO2 當 Mask**："
                       f"在此製程/溫度條件下，SiO2 被蝕刻的速度比氮化矽慢，"
                       f"所需 Mask 厚度較薄（{sio2_needed_nm:.0f} nm vs. "
                       f"{nitride_needed_nm:.0f} nm）。")
            else:
                st.info("兩種 Mask 材料在此條件下所需厚度相近，可依製程相容性"
                       "（如應力、沉積方式）另行選擇。")

            if gc_row["pplus_etch_stop"] == "Yes":
                st.caption("附註：此製程對 p++ 重摻雜矽有蝕刻終止（etch-stop）效果，"
                          "若結構設計允許，也可考慮用 p++ 摻雜層取代或輔助 Mask。")
        else:
            st.warning("蝕刻速率為 0，無法計算所需時間。")

    else:  # 給定 Mask 厚度反算最大深度
        mask_material = st.selectbox("選擇 Mask 材料", ["SiO2", "氮化矽 (Si3N4)"])
        mask_thickness_nm = st.number_input("現有 Mask 厚度 (nm)", min_value=0.0,
                                             value=500.0, step=50.0)
        mrate = sio2_rate_nm_min if mask_material == "SiO2" else nitride_rate_nm_min

        if mrate > 0:
            max_time_min = mask_thickness_nm / (mrate * safety_factor)
            max_depth_um6 = si_rate_um_min * max_time_min
            st.success(f"此 Mask 厚度（含安全係數）最多可承受 "
                       f"**{max_time_min:.1f} 分鐘**的蝕刻，"
                       f"對應最大 Si 蝕刻深度約 **{max_depth_um6:.1f} µm**。")
        else:
            st.info(f"{mask_material} 在此製程條件下幾乎不被蝕刻"
                    "（參考速率 ≈ 0），理論上蝕刻時間不受此 Mask 限制。")

    st.markdown("---")
    st.subheader("各製程 Mask 材料操作點蝕刻速率一覽")
    st.dataframe(build_general_char_display(), use_container_width=True)

# =========================================================
# ⑦ 參考文獻與資料來源
# =========================================================
elif page == "⑦ 參考文獻與資料來源":
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
- Tan et al. -- 引用自二次文獻（見下方 Handbook），
  原始期刊出處尚待進一步查證。
- ~~Wind et al. (37 wt%)~~ -- 已於本版移除。原資料組在測試中被發現
  數值異常（{100}/{110}/{111} 三個晶面在常用溫度下算出的選擇比明顯
  偏離其他來源與講義數值），且其原始期刊出處始終無法查證，故直接
  從資料庫移除，不再提供使用者選擇，避免誤用。

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

### 製程操作點資料（④、⑥模組使用）
- 課堂 PDF《Etching》p.2 "General Characteristics of Silicon Etching
  Operations"。原始出處：Kovacs, G.T.A., Maluf, N.I., Peterson, K.E.
  (1998). "Bulk micromachining of silicon." *Proceedings of the IEEE*,
  86(8), 1536-1551.
- 課堂 PDF《Etching》p.7 "Comparison of Etch Rates for Selected
  Etchants and Target Materials"。原始出處：Williams, K. & Muller, R.
  （作者全名、發表年份見 PDF 原文，本應用程式僅轉引課堂投影片內容，
  建議查證原始論文以取得完整書目資訊）。

### 異向性蝕刻幾何（⑤模組使用）
- 課堂 PDF《Etching》"Anisotropic Etching" 頁示意圖：{111} 面與晶圓
  表面夾角 54.7°，兩種異向性蝕刻方式（orientation-dependent etching,
  ODE；vertical etching）。54.7° 為矽晶格中 {111} 與 {100} 晶面夾角的
  標準幾何值（arccos(1/√3) ≈ 54.74°）。

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
2. **Tan et al. / Tabata et al.** 的原始期刊出處：
   目前僅能確認其數據被 Gosálvez et al. (2015) 的 Table 22.3 引用，
   原始期刊、卷期、頁碼尚未逐一查證。（Wind et al. 37 wt% 資料組已於
   測試中發現異常並移除，見上方 KOH 動力學參數說明。）
3. **EDP 的前置因子 R0**：文獻 (Dutta et al. 2011) 僅報告活化能
   0.34 eV，未提供 R0 數值，本程式的 R0 是用課堂 PDF 單點資料反推。
4. **HF:HNO3:CH3COOH（HNA系統）**：尚未查到公開發表的 Arrhenius
   參數 (Ea, R0)，目前完全未納入動力學計算器，僅在④、⑥模組以 PDF
   表格中的單一操作點數值（25°C, 1-20 µm/min）呈現。
5. **PDF 講義的「100:1」選擇比 vs. 計算器算出的比值**：講義 General
   Characteristics 表格給的 {111}/{100}=100:1 是**橫跨多種操作條件的
   典型數量級概估值**，而②計算器套用的是 Seidel et al. (1990) 在
   **34 wt%、特定溫度**下的精確 Arrhenius 方程式，兩者本來就不是同一
   件事：KOH 的 {111}/{100} 選擇比對濃度非常敏感，濃度越高選擇比通常
   越低（約 20 wt% 時可達 100:1 以上，34 wt% 高濃度下常落在 50:1～
   70:1 附近），所以計算器在 34 wt%、80°C 算出約 60:1 是合理的實驗
   範圍，跟講義「100:1」沒有矛盾，只是後者是概估值、前者是精確模型
   在特定條件下的結果。詳細說明已加在②頁面的展開框中。
5. **④、⑥模組的材料選擇比資料（PDF p.7 表格）**：為課堂投影片轉引，
   原始論文完整書目資訊（作者全名、期刊、年份）待查證，含 "<"/">"
   的數值以估計值處理，使用時請對照原始 PDF 或原始文獻確認。
""")
