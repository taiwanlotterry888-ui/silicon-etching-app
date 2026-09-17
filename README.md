# Silicon Etching Process Evaluator（矽蝕刻製程評估系統）

依據課堂 PDF《Etching》內容製作的網頁應用，協助評估矽蝕刻製程的蝕刻速率、
材料選擇比與蝕刻剖面。使用 **Streamlit** 框架，純 Python 撰寫。

## 檔案結構
```
.
├── app.py              # 主程式（所有功能模組都在這裡）
├── requirements.txt     # 相依套件
└── README.md
```

## 功能模組
1. 濕蝕刻速率計算器（KOH / TMAH / EDP / HF:HNO3:CH3COOH）
2. 蝕刻劑材料選擇比查詢（依蝕刻劑查/依材料反查）
3. 異向性蝕刻剖面模擬器（{111}面 54.7° 幾何繪圖）
4. 乾蝕刻（電漿）計算器（SF6 / CF4+CHF3+He / O2 plasma）
5. 蝕刻時間反算 / Mask 材料建議工具

## 本機測試方式
```bash
pip install -r requirements.txt
streamlit run app.py
```
瀏覽器會自動開啟 `http://localhost:8501`

## 部署到 GitHub + 產生網址（Streamlit Community Cloud，免費）

### 步驟一：上傳到 GitHub
1. 到 GitHub 建立一個新的 repository（例如 `silicon-etching-app`）
2. 把 `app.py`、`requirements.txt`、`README.md` 三個檔案上傳（或用 git push）：
   ```bash
   git init
   git add app.py requirements.txt README.md
   git commit -m "Initial commit: silicon etching evaluator"
   git branch -M main
   git remote add origin https://github.com/<你的帳號>/silicon-etching-app.git
   git push -u origin main
   ```

### 步驟二：部署到 Streamlit Community Cloud
1. 前往 https://share.streamlit.io/ 並用 GitHub 帳號登入
2. 點選 **"New app"**
3. 選擇剛剛的 repository、branch（main）
4. **Main file path** 填寫 `app.py`
5. 點選 **Deploy**，約 1~2 分鐘後會產生一個網址，格式類似：
   ```
   https://<你的帳號>-silicon-etching-app.streamlit.app
   ```
6. 之後只要 `git push` 更新程式碼，Streamlit Cloud 會自動重新部署。

## 資料來源
所有蝕刻速率、選擇比數據皆整理自課堂 PDF《Etching》：
- General Characteristics of Silicon Etching Operations
  （Kovacs, Maluf & Peterson, "Bulk micromachining of silicon",
  Proceedings of the IEEE, v.86, No.8, 1998）
- Comparison of Etch Rates for Selected Etchants and Target Materials
  （K. Williams & R. Muller）

> 實際製程中蝕刻速率會受溫度、濃度、曝露面積與微結構影響而有所差異，
> 本工具僅供教學與初步製程評估參考。
