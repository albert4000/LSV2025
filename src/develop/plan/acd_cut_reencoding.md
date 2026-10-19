# 計劃：Pipeline cut 的 ACD 重編碼——可行性量測、實作與驗證

## 0. 給 Claude Code 的說明

這份計劃要回答一個問題：**在有延遲限制的 min-register cut 上，用 Ashenhurst-Curtis decomposition（ACD）式的重編碼，能不能降低總面積？**

工作分成四個階段，**Phase 1 結束後有 go/no-go 判斷**，沒過門檻的電路不進 Phase 2。請依序做，每個階段結束時停下來回報結果再繼續。

開始前請先讀現有的 min-cut 程式，確認下列假設；任何一項不符就先回報，不要自行改假設。

### 假設（需確認）

| 項目 | 假設 |
|---|---|
| 表示法 | AIG，benchmark 為 EPFL combinational suite（arithmetic + random/control） |
| Delay model | unit delay，一個 AND node = 1 level |
| Stage 數 | 2 級（單一 cut），T = 每級允許的最大 level 數 |
| 面積模型 | 1 個 AND node ≈ 1 個 NAND2；1 個 register = 10 個 NAND2 |
| 目標函數 | 總面積 = AND node 數 + 10 × register 數，在兩級深度都 ≤ T 的條件下最小化 |
| 現有工具 | 已有 delay-constrained min-cut，可輸出 cut 上的 node 集合 |

### 名詞

- **cut 訊號集合 C**：輸出要被打 register 的 node（含直接穿越的 PI）。
- **group G ⊆ C**：一起被重編碼的一群 cut 訊號。
- **window leaves L**：G 上游的一組 node，G 的每個訊號都是 L 的函數，|L| ≤ K。
- **image Img(G)**：L 取遍所有值時，G 實際出現的相異向量集合。
- **r(G) = ⌈log₂|Img(G)|⌉**：G 重編碼後所需的最少線數。
- **upstream slack(G)** = T − max(G 中訊號的 arrival level)。
- **downstream slack(G)** = T − max(G 中訊號到 PO 的最長 level 數)。

## 1. 背景與待驗證的假說

現有量測顯示，許多電路在建議點的 cut 寬度已超過 PI 數（例如 log2 199 vs 32、sin 66 vs 24、arbiter 513 vs 256），代表 cut 上的訊號有冗餘。但「把 cut 壓回 PI 數」永遠有一個 trivial 解：把 register 移到輸入端，而這會違反延遲限制。所以真正的問題是 **在 slack 內能否算出夠便宜的編碼**。

以下是待驗證的假說，**不是已知事實**，請用資料確認或推翻：

- **H1**：控制型、one-hot 輸出的電路（arbiter、priority、router、voter）有大的可壓縮量，且 encoder/decoder 便宜。
- **H2**：算術 datapath（hyp、sqrt、sin、log2、multiplier、div、square）雖然 cut > PI 數，但冗餘來自 carry-save 之類的表示法，壓縮所需的深度超過 slack，淨收益接近零。
- **H3**：max 的 cut 是兩個候選值，幾乎不可壓縮。
- **H4**：對有膝點的電路，把 T 從 50% 放寬到建議值所省的 register，比 ACD 在 T=50% 能省的多。

## 2. Phase 0：重現 baseline

1. 讀懂現有 min-cut 程式的介面、delay model、輸出格式。
2. 對下表 20 個電路、T ∈ {50%, 建議 T} 重跑，確認 register 數與下表一致。
3. 對每個 (電路, T) 輸出 cut 訊號清單，附每個訊號的 arrival level 與到 PO 的最長 level 數。

| 電路 | 深度 | reg @ T=50% | 建議 T | reg @ 建議 T | 不加限制 |
|---|---|---|---|---|---|
| bar | 12 | 259 | 7 | 131 | 128 |
| hyp | 24801 | 306 | 13641 | 257 | 128 |
| sqrt | 5058 | 166 | 2782 | 137 | 64 |
| sin | 225 | 86 | 124 | 66 | 24 |
| mem_ctrl | 114 | 1554 | 63 | 1444 | 1193 |
| multiplier | 274 | 221 | 151 | 190 | 128 |
| log2 | 444 | 210 | 245 | 199 | 32 |
| arbiter | 87 | 513 | 44 | 513 | 129 |
| max | 287 | 258 | 144 | 258 | 130 |
| dec | 3 | 16 | 2 | 16 | 8 |
| adder | 255 | 193 | 192 | 161 | 129 |
| square | 250 | 192 | 188 | 159 | 64 |
| div | 4372 | 335 | 3279 | 299 | 128 |
| priority | 250 | 72 | 188 | 40 | 8 |
| router | 54 | 50 | 41 | 22 | 3 |
| voter | 70 | 86 | 53 | 16 | 1 |
| cavlc | 16 | 63 | 12 | 28 | 10 |
| int2float | 16 | 26 | 12 | 13 | 7 |
| ctrl | 10 | 27 | 8 | 10 | 7 |
| i2c | 20 | 169 | 15 | 145 | 141 |

**完成條件**：register 數全部吻合；不吻合的項目列出差異與原因。

## 3. Phase 1：可壓縮量量測（不修改電路）

目的是用便宜的方法估計「最多能省幾個 register」，決定哪些電路值得進 Phase 2。

### 3.1 全域取樣（粗略指標）

- 對每個 (電路, T)：PI 數 ≤ 20 時窮舉，否則模擬 2²⁰ 個隨機 pattern（bit-parallel）。
- 用 hash 統計 cut 向量的相異數 N，回報 ⌈log₂N⌉。
- **注意**：取樣只能給 N 的下界。若 N 等於取樣數，標記為 `saturated`，此指標對該電路無效，以 3.2 為準。

### 3.2 Windowed image 分析（主要指標）

1. **分群**：把 C 切成互斥的 group。貪婪法：取一個未分配的 cut 訊號當種子，反覆加入與目前 group 的 window leaves 重疊最多的 cut 訊號，直到 leaves 數會超過 K。
2. **取 window**：對每個 group 用 reconvergence-driven cut 往上游找 leaves L，|L| ≤ K，預設 K = 16。leaves 不必是 PI。
3. **算 image**：對 L 窮舉 2^|L| 個值，統計 G 的相異向量數，得 r(G)。把 L 當成自由變數會高估 image，所以結果是保守的（不會謊報可壓縮量）。
4. **排除 trivial 解**：若 r(G) 的壓縮只能靠「把 register 移到 L」達成，且 L 的 arrival 讓下游深度超過 T，這個節省在 Phase 2 會被 slack 檢查擋掉。Phase 1 先照實記錄，另外標記 `r(G) == |L|` 的 group。
5. **函數相依檢查**：找 G 的最小子集 S，使 G 其餘訊號都是 S 的函數（貪婪刪除即可）。這種情況不需要 encoder，只需下游 decoder，是最便宜的壓縮，單獨統計。
6. 對 K ∈ {8, 12, 16} 各跑一次，觀察節省量對 window 大小的敏感度。

### 3.3 輸出

`phase1_groups.csv`，每個 group 一列：

```
circuit, T, group_id, |G|, |L|, |Img|, r, saving(=|G|-r),
fd_subset_size, upstream_slack, downstream_slack, r_equals_L
```

`phase1_summary.csv`，每個 (電路, T) 一列：

```
circuit, T, cut_size, sum_saving, sum_saving_fd_only,
sum_saving_with_slack>0, global_log2N, saturated
```

### 3.4 Go/no-go 門檻

某電路進 Phase 2 的條件（兩者皆須成立）：

- `sum_saving_with_slack>0` ≥ cut_size 的 5%，且 ≥ 4 個 register；
- 至少一個 group 的 upstream 或 downstream slack ≥ 2 level。

**完成條件**：產出兩個 CSV，加一段文字對 H1–H3 逐條說明資料支持或不支持。在此停下回報。

## 4. Phase 2：合成 encoder / decoder

只對通過門檻的電路做。對每個 group 依下列順序嘗試，取第一個被接受的方案：

1. **函數相依**：只 register S，下游用 decoder 重建 G \ S。
2. **一般重編碼**：
   - 為 Img(G) 的每個值指定 r-bit code。先試「盡量沿用 G 中既有訊號當 code bit」，不足的 bit 再指定新 code。
   - Encoder：L → r bits，由 truth table 合成。
   - Decoder：r bits → G，未使用的 code 當 don't care。
   - 用 ABC 合成並最佳化（truth table → `strash` → `resyn2` 或 `dc2`），取最小 AND 數的結果。

### 接受條件（全部成立才套用）

- 面積：10 × (|G| − r) > AND(encoder) + AND(decoder) − 因此變成 dangling 而移除的 AND 數。
- 上游延遲：encoder 輸出的 arrival level ≤ T。
- 下游延遲：decoder 插入後，第二級深度 ≤ T。
- 套用後重新計算全域 level，確認沒有其他路徑因共用邏輯而變長。

group 之間可能共用上游邏輯，請逐一套用並在每次套用後更新 level 與 AND 數，不要一次批次估算。

### 輸出

`phase2_results.csv`：

```
circuit, T, reg_before, reg_after, and_before, and_after,
area_before, area_after, area_delta_pct,
groups_tried, groups_accepted, groups_rejected_area, groups_rejected_delay
```

## 5. Phase 3：驗證

每個修改過的電路都必須通過：

1. **功能等價**：把 register 當成導線，與原始 AIG 做 `cec`。不等價即視為失敗，不得回報面積數字。
2. **延遲**：兩級深度各自 ≤ T，獨立重新計算，不沿用 Phase 2 的增量結果。
3. **register 數**：從最終 netlist 重新數，與 Phase 2 回報值一致。
4. **面積**：重新計算 AND 數 + 10 × register 數。

### 單元測試（先寫，再實作 Phase 1/2）

| 測試 | 預期 |
|---|---|
| 手寫 3-to-8 one-hot decoder，cut 在 one-hot 輸出 | image = 8，r = 3 |
| 手寫電路：cut 上有訊號 a、b、a∧b | 函數相依檢查找出 S = {a, b} |
| cut 訊號互相獨立（直接穿越的 PI） | saving = 0 |
| dec @ T=2 | Phase 1 回報可壓到 8；Phase 2 因下游 slack 不足而全數拒絕 |
| 刻意塞入錯誤 decoder | `cec` 報不等價，流程中止 |

## 6. Phase 4：與放寬 T 比較（H4）

對每個電路整理三個數字，放在同一張表：

- A：T=50%，原始 min-cut 的總面積
- B：T=50%，加 ACD 重編碼後的總面積
- C：T=建議值，原始 min-cut 的總面積

回報 (A−B)/A 與 (A−C)/A。若建議 T 下也跑了 ACD，再加一欄 D。

## 7. 最終交付

- 上述 CSV 檔與產生它們的 script，可一鍵重跑。
- `REPORT.md`：
  - H1–H4 逐條結論，附支持的數字；
  - 哪些電路有淨面積改善、改善多少；
  - 被拒絕的 group 中，面積不划算與延遲不夠各佔多少；
  - 沒做到或結果存疑的部分，明確列出。
- 所有修改後的 AIG 檔，以及對應的 `cec` log。

## 8. 範圍外（這次不做）

- 考慮下游 observability don't care 的完整 column multiplicity（比 image 更強，但需要下游窮舉）。若 Phase 1 顯示 image 法在某些電路差一點就划算，再另外評估。
- 多於 2 級的 pipeline。
- 架構層級的重寫（例如把 ripple adder 改成 carry-select）。
- Technology mapping 後的真實面積與時序。

## 9. 執行注意事項

- hyp（深度 24801）、div、sqrt 的電路很大，Phase 1 請先在小電路（cavlc、int2float、ctrl、dec、router）上跑通再擴大。
- 單一電路超過 30 分鐘就中止並回報，不要默默縮小 K 或取樣數。
- 任何與本計劃假設不符的發現，優先回報，不要為了讓流程跑完而繞過。