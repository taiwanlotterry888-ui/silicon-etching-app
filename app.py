# -*- coding: utf-8 -*-
"""
Silicon Etching Process Evaluator
-----------------------------------
一个用于评估矽（Silicon）蝕刻製程的網頁應用。
資料來源：課堂 PDF「Etching」中的兩張表格
 (1) General Characteristics of Silicon Etching Operations
 (2) Comparison of Etch Rates for Selected Etchants and Target Materials
     (Source: Kovacs, Maluf & Peterson, IEEE Proceedings, 1998;
      Williams & Muller)

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

# =========================================================
# 1. 內建資料庫（依 PDF 表格整理）
# =========================================================

# ---- Table 1: General Characteristics of Silicon Etching Operations ----
general_wet = pd.DataFrame([
    {"Etchant": "HF:HNO3:CH3COOH", "Category": "Wet",
     "Temp_C_min": 25, "Temp_C_max": 25,
     "EtchRate_um_min_min": 1, "EtchRate_um_min_max": 20,
     "Selectivity_111_100": None,
     "Nitride_nm_min": None, "Nitride_note": "Low",
     "SiO2_nm_min": 10, "SiO2_nm_max": 30,
     "p++_stop": "No"},
    {"Etchant": "KOH", "Category": "Wet",
     "Temp_C_min": 70, "Temp_C_max": 90,
     "EtchRate_um_min_min": 0.5, "EtchRate_um_min_max": 2,
     "Selectivity_111_100": "100:1",
     "Nitride_nm_min": 1, "Nitride_note": "<1",
     "SiO2_nm_min": 10, "SiO2_nm_max": 10,
     "p++_stop": "Yes"},
    {"Etchant": "Ethylene-diamine pyrocatechol (EDP)", "Category": "Wet",
     "Temp_C_min": 115, "Temp_C_max": 115,
     "EtchRate_um_min_min": 0.75, "EtchRate_um_min_max": 0.75,
     "Selectivity_111_100": "35:1",
     "Nitride_nm_min": 0.1, "Nitride_note": "0.1",
     "SiO2_nm_min": 0.2, "SiO2_nm_max": 0.2,
     "p++_stop": "Yes"},
    {"Etchant": "N(CH3)4OH (TMAH)", "Category": "Wet",
     "Temp_C_min": 90, "Temp_C_max": 90,
     "EtchRate_um_min_min": 0.5, "EtchRate_um_min_max": 1.5,
     "Selectivity_111_100": "50:1",
     "Nitride_nm_min": 0.1, "Nitride_note": "<0.1",
     "SiO2_nm_min": 0.1, "SiO2_nm_max": 0.1,
     "p++_stop": "Yes"},
    {"Etchant": "SF6 (Dry/Plasma)", "Category": "Dry",
     "Temp_C_min": 0, "Temp_C_max": 100,
     "EtchRate_um_min_min": 0.1, "EtchRate_um_min_max": 0.5,
     "Selectivity_111_100": None,
     "Nitride_nm_min": 200, "Nitride_note": "200",
     "SiO2_nm_min": 10, "SiO2_nm_max": 10,
     "p++_stop": "No"},
    {"Etchant": "SF6/C4F8 (DRIE)", "Category": "Dry",
     "Temp_C_min": 20, "Temp_C_max": 80,
     "EtchRate_um_min_min": 1, "EtchRate_um_min_max": 3,
     "Selectivity_111_100": None,
     "Nitride_nm_min": 200, "Nitride_note": "200",
     "SiO2_nm_min": 10, "SiO2_nm_max": 10,
     "p++_stop": "No"},
])

# ---- Table 2: Isotropic / Dry Etchants vs Target material etch rate (nm/min) ----
# 欄位: 材料的蝕刻速率 (nm/min)。 None 代表 PDF 原表中為 "-"（無數據）。
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

# {111} face angle for <100> silicon wafers (constant used in KOH/TMAH/EDP anisotropic etching)
ANGLE_111 = 54.7  # degrees


# =========================================================
# Sidebar 導覽
# =========================================================
st.sidebar.title("🔬 Silicon Etching Evaluator")
page = st.sidebar.radio(
    "選擇功能模組",
    [
        "① 首頁說明",
        "② 濕蝕刻速率計算器",
        "③ 蝕刻劑材料選擇比查詢",
        "④ 異向性蝕刻剖面模擬器",
        "⑤ 乾蝕刻（電漿）計算器",
        "⑥ 蝕刻時間反算 / Mask建議",
    ],
)

# =========================================================
# ① 首頁
# =========================================================
if page == "① 首頁說明":
    st.title("矽蝕刻製程評估系統")
    st.markdown("""
本應用依據課堂 PDF「Etching」內容建置，目的是協助工程師快速評估矽蝕刻製程的
**蝕刻速率、材料選擇比與蝕刻輪廓**，涵蓋濕蝕刻（Wet Etching）與乾蝕刻
（Dry/Plasma Etching）兩大類。

### 功能模組
| 模組 | 功能 | 對應 PDF 內容 |
|---|---|---|
| 濕蝕刻速率計算器 | 依蝕刻劑/溫度查詢蝕刻速率、{111}/{100}選擇比、p++ etch stop | General Characteristics 表 |
| 蝕刻劑材料選擇比查詢 | 查詢特定蝕刻劑對各材料(SiO2/Nitride/PSG/Al/Ti/PR)的蝕刻速率並算選擇比 | Isotropic Etchants 表 |
| 異向性蝕刻剖面模擬器 | 依 Mask 開口寬度、蝕刻深度模擬 KOH/TMAH 異向性蝕刻的剖面（54.7°） | Anisotropic Etching 頁 |
| 乾蝕刻計算器 | 依氣體/功率查詢蝕刻速率 | Dry Etchants 表 |
| 蝕刻時間反算/Mask建議 | 輸入目標深度反算所需時間；依選擇比建議 Mask 材料 | 綜合應用 |

> ⚠️ 資料皆為 PDF 課程投影片所附之參考值，實際蝕刻速率會因設備、濃度、溫度、
> 曝露面積與微結構而異（如 PDF 表格備註所述）。
""")

# =========================================================
# ② 濕蝕刻速率計算器
# =========================================================
elif page == "② 濕蝕刻速率計算器":
    st.header("濕蝕刻 / 乾蝕刻 速率計算器")
    st.caption("資料來源：General Characteristics of Silicon Etching Operations")

    etchant_choice = st.selectbox("選擇蝕刻劑", general_wet["Etchant"].tolist())
    row = general_wet[general_wet["Etchant"] == etchant_choice].iloc[0]

    col1, col2 = st.columns(2)
    with col1:
        st.metric("蝕刻類別", row["Category"])
        st.metric("建議溫度範圍 (°C)",
                  f"{row['Temp_C_min']} ~ {row['Temp_C_max']}")
        st.metric("蝕刻速率範圍 (µm/min)",
                  f"{row['EtchRate_um_min_min']} ~ {row['EtchRate_um_min_max']}")
    with col2:
        st.metric("{111}/{100} 選擇比",
                  row["Selectivity_111_100"] if row["Selectivity_111_100"] else "無資料")
        st.metric("p++ Etch Stop", row["p++_stop"])
        st.metric("SiO2 蝕刻速率 (nm/min)",
                  f"{row['SiO2_nm_min']} ~ {row['SiO2_nm_max']}")

    st.markdown("---")
    st.subheader("目標深度 → 預估時間")
    target_depth_um = st.number_input("目標蝕刻深度 (µm)", min_value=0.0,
                                       value=100.0, step=10.0)
    rate_choice = st.slider(
        "使用蝕刻速率 (µm/min)（可在範圍內調整以模擬濃度/溫度差異）",
        float(row["EtchRate_um_min_min"]), float(row["EtchRate_um_min_max"]),
        float(row["EtchRate_um_min_min"]))
    if rate_choice > 0:
        est_time_min = target_depth_um / rate_choice
        st.success(f"預估所需時間：約 **{est_time_min:.1f} 分鐘** "
                    f"（{est_time_min/60:.2f} 小時）")

    st.markdown("---")
    st.subheader("完整資料表")
    st.dataframe(general_wet, use_container_width=True)

# =========================================================
# ③ 蝕刻劑材料選擇比查詢
# =========================================================
elif page == "③ 蝕刻劑材料選擇比查詢":
    st.header("蝕刻劑 × 材料 蝕刻速率／選擇比查詢")
    st.caption("資料來源：Comparison of Etch Rates for Selected Etchants and Target Materials")

    tab1, tab2 = st.tabs(["依蝕刻劑查詢", "依材料反查（找出蝕刻該材料最快/最慢的蝕刻劑）"])

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
        sub = etchant_material_rate[["Etchant", "Type", material_choice]].copy()
        sub = sub.rename(columns={material_choice: "EtchRate_nm_min"})
        sub = sub.dropna(subset=["EtchRate_nm_min"]).sort_values(
            "EtchRate_nm_min", ascending=False)
        st.write(f"針对 **{material_choice}**，各蝕刻劑速率（由快到慢）：")
        st.dataframe(sub, use_container_width=True)
        if len(sub) > 0:
            st.success(f"蝕刻速度最快：**{sub.iloc[0]['Etchant']}** "
                        f"({sub.iloc[0]['EtchRate_nm_min']} nm/min)")
            st.info(f"蝕刻速度最慢（最適合當作對此材料無害的 Mask）："
                    f"**{sub.iloc[-1]['Etchant']}** "
                    f"({sub.iloc[-1]['EtchRate_nm_min']} nm/min)")

    st.markdown("---")
    st.subheader("完整資料表")
    st.dataframe(etchant_material_rate, use_container_width=True)

# =========================================================
# ④ 異向性蝕刻剖面模擬器
# =========================================================
elif page == "④ 異向性蝕刻剖面模擬器":
    st.header("異向性蝕刻（Anisotropic Etching）剖面模擬")
    st.caption("模擬 <100> 矽晶片以 KOH/TMAH/EDP 進行 orientation-dependent etching (ODE) 的剖面")

    col1, col2, col3 = st.columns(3)
    with col1:
        opening_um = st.number_input("Mask 開口寬度 (µm)", min_value=1.0,
                                      value=100.0, step=5.0)
    with col2:
        depth_um = st.number_input("蝕刻深度 (µm)", min_value=1.0,
                                    value=50.0, step=5.0)
    with col3:
        wafer_thickness_um = st.number_input("晶片厚度 (µm，選填，用於判斷是否貫穿)",
                                              min_value=0.0, value=0.0, step=10.0)

    mode = st.radio("蝕刻模式", ["異向性 (ODE, 54.7° {111}面)", "等向性 (Isotropic)"])

    # 幾何計算：異向性蝕刻時，底部寬度隨深度縮減
    # 側壁與垂直方向夾角 = 90 - 54.7 = 35.3°；每邊縮減量 = depth / tan(54.7°)
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
        # 等向性：假設 undercut ≈ 蝕刻深度（等速蝕刻各方向）
        shrink_per_side = depth_um
        bottom_width = opening_um  # 底部維持開口寬（但有underercut）
        v_groove = False
        st.info(f"等向性蝕刻：Mask 下方 undercut ≈ **{shrink_per_side:.1f} µm**（各方向蝕刻速率相同）")

    if wafer_thickness_um > 0:
        if depth_um >= wafer_thickness_um:
            st.warning("⚠️ 蝕刻深度已達到或超過晶片厚度，可能已貫穿晶片！")

    # 繪圖
    fig, ax = plt.subplots(figsize=(6, 4))
    top_y = 0
    bottom_y = -depth_um
    half_open = opening_um / 2

    if mode.startswith("異向性"):
        if v_groove:
            apex_depth = half_open * np.tan(np.radians(ANGLE_111))
            xs = [-half_open, 0, half_open]
            ys = [top_y, -apex_depth, top_y]
        else:
            half_bottom = bottom_width / 2
            xs = [-half_open, -half_bottom, half_bottom, half_open]
            ys = [top_y, bottom_y, bottom_y, top_y]
        ax.plot(xs, ys, 'b-', linewidth=2)
        ax.fill(xs, ys, alpha=0.3, color='deepskyblue')
    else:
        half_bottom = half_open + shrink_per_side  # undercut 使底部變寬（示意）
        xs = [-half_open, -half_bottom, half_bottom, half_open]
        ys = [top_y, bottom_y, bottom_y, top_y]
        ax.plot(xs, ys, 'g-', linewidth=2)
        ax.fill(xs, ys, alpha=0.3, color='lightgreen')

    # Mask 示意
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

    st.markdown("""
**計算原理：**
- 異向性蝕刻：矽 <100> 晶片沿 {111} 面蝕刻，側壁角度固定為 **54.7°**，
  每側內縮量 = 深度 ÷ tan(54.7°)。
- 等向性蝕刻：假設各方向蝕刻速率相同，Mask 下方會產生 undercut，
  undercut 量 ≈ 蝕刻深度。
""")

# =========================================================
# ⑤ 乾蝕刻（電漿）計算器
# =========================================================
elif page == "⑤ 乾蝕刻（電漿）計算器":
    st.header("乾蝕刻 / 電漿蝕刻 速率計算器")
    st.caption("資料來源：Dry Etchants (CF4+CHF3+He, SF6+He, SF6, O2 plasma)")

    dry_df = etchant_material_rate[etchant_material_rate["Type"] == "Dry"]
    dry_choice = st.selectbox("選擇乾蝕刻條件（氣體/功率）",
                               dry_df["Etchant"].tolist())
    row = dry_df[dry_df["Etchant"] == dry_choice].iloc[0]

    material_choice = st.selectbox("選擇欲蝕刻材料", MATERIAL_COLS, index=0)
    rate = row[material_choice]

    if pd.isna(rate):
        st.warning("此條件對所選材料無資料。")
    else:
        st.metric(f"{material_choice} 蝕刻速率", f"{rate} nm/min")

        target_nm = st.number_input("目標蝕刻深度 (nm)", min_value=0.0,
                                     value=500.0, step=50.0)
        if rate > 0:
            est_time = target_nm / rate
            st.success(f"預估所需時間：約 **{est_time:.1f} 分鐘**")
        else:
            st.info("此條件下該材料幾乎不被蝕刻（速率 = 0），可作為抗蝕刻層。")

    st.markdown("---")
    st.subheader("完整乾蝕刻資料表")
    st.dataframe(dry_df, use_container_width=True)

# =========================================================
# ⑥ 蝕刻時間反算 / Mask 建議
# =========================================================
elif page == "⑥ 蝕刻時間反算 / Mask建議":
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

        target_material = row["TargetMaterial"]
        target_rate = None
        # 嘗試找出與 TargetMaterial 名稱對應的欄位速率（若在 MATERIAL_COLS 中）
        for m in MATERIAL_COLS:
            if m.lower().replace("_", " ") in target_material.lower() or \
               target_material.lower() in m.lower():
                target_rate = row[m]
                break

        if best_mask[1] > 0 and target_rate and target_rate > 0:
            selectivity = target_rate / best_mask[1]
            st.info(f"與目標蝕刻材料的選擇比 ≈ **{selectivity:.1f} : 1**")
        elif best_mask[1] == 0:
            st.info("此材料幾乎完全不被蝕刻 → 選擇比視為極高，是理想的 Mask 材料。")

    st.markdown("---")
    st.subheader("跨模組小工具：目標深度 → 建議蝕刻劑（速率排序）")
    material_target = st.selectbox("選擇欲蝕刻的材料", MATERIAL_COLS, key="m2")
    target_depth_nm = st.number_input("目標蝕刻深度 (nm)", min_value=0.0,
                                       value=1000.0, step=100.0, key="d2")

    sub = etchant_material_rate[["Etchant", "Type", material_target]].dropna()
    sub = sub.rename(columns={material_target: "Rate"})
    sub = sub[sub["Rate"] > 0].copy()
    if len(sub) > 0:
        sub["EstTime_min"] = target_depth_nm / sub["Rate"]
        sub = sub.sort_values("EstTime_min")
        st.dataframe(sub, use_container_width=True)
        st.success(f"最快方案：**{sub.iloc[0]['Etchant']}**，"
                    f"預估 {sub.iloc[0]['EstTime_min']:.1f} 分鐘")
    else:
        st.warning("查無可蝕刻此材料且速率>0的資料。")

st.sidebar.markdown("---")
st.sidebar.caption("資料來源：課堂 PDF《Etching》\n"
                    "Kovacs, Maluf & Peterson (1998); Williams & Muller")
