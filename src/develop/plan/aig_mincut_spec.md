# AIG Min-Cut：固定 stage 數的 pipeline register 最小化（實作規格）

## 1. 目標

輸入一個組合 AIG 與 stage 數 S，在 AIG 上插入 S−1 層 register cut，使得：

1. **Balanced**：每條 PI→PO 路徑恰好穿過 S−1 個 register。
2. **最小化 register 總數**（register 面積 = 10 × NAND2，見第 6 節的權重）。
3. （可選）每個 stage 的組合延遲 ≤ T。

本階段**不做**邏輯複製、encoder/decoder、Curtis 重編碼。只做結構性 min-cut。輸出的 cut 要能給後續階段使用。

## 2. 記號與模型

- AIG：節點集合 V（PI、AND、可有 constant），邊 u→v 表示 u 是 v 的 fanin。反相邊（complemented edge）**不影響 cut**，register 存節點的正相值，反相是免費的。
- PO：每個 PO 有一個 driver 節點。
- 每個節點 v 指派 stage 編號 s(v) ∈ {0..S−1}，PI 的 s = 0，PO driver 需要被 register 到 stage S−1。
- 一條邊 u→v 跨過 s(v) − s(u) 個 cut。
- **Register 共用**：同一個節點 u 的所有 fanout 共用一條 register chain，長度為
  `m(u) − s(u)`，其中 `m(u) = max( s(w) for w in fanout(u) ∪ {PO sink if u 是 PO driver, 值 S−1} )`。
- **成本 = Σ_u W(u) · (m(u) − s(u))**，W(u) 為該節點一個 register 的面積權重（預設 1，面積模式下為 10）。
  注意是對「節點」計算，不是對邊計算。

## 3. 演算法總覽

| 階段 | 問題 | 方法 | 精確性 |
|---|---|---|---|
| A | S=2（單一 cut） | node splitting + max-flow / min-cut | 精確 |
| B | 一般 S | s(v)、m(u) 的差分約束 LP，對偶為 min-cost flow | 精確（TU，LP 解為整數） |
| C | 時序限制 T | level 窗口 + lazy constraint generation | 啟發式修正 |
| D | 輸出 | cut 集合、報表、驗證 | — |

建議實作順序：**A → 驗證 → B → 驗證 → C**。A 是 B 在 S=2 的特例，可以用來交叉驗證 B。

## 4. 階段 A：單一 cut（S = 2）

### 4.1 建圖

對每個 AIG 節點 v 拆成 `v_in`、`v_out`：

| 邊 | 容量 | 意義 |
|---|---|---|
| `v_in → v_out` | W(v) | 在 v 的輸出放 register（v 被 cut）。這條邊被切斷 = 付一份 register 成本 |
| `u_out → v_in`（對每條 AIG 邊 u→v） | ∞ | 結構連線 |
| `v_in → u_in`（對每條 AIG 邊 u→v） | ∞ | **反向邊，強制 cut 為 closed**：若 v 在 stage 0，則所有 fanin 也必須在 stage 0，避免路徑穿過 cut 兩次 |
| `source → PI_in` | ∞ | |
| `PO_driver_out → sink` | ∞ | PO driver 必須被 cut 或位於 cut 之後 |

- Constant 節點：不連入圖（不需要 register）。
- 若 PO driver 就是 PI，PI 本身會被 cut，表示在輸入端放 register，這是合法解。

### 4.2 求解

- 跑 max-flow（Dinic，大圖用 push-relabel 或 BK），min-cut 值 = 最少 register 數（W=1 時）。
- 從 residual graph 取 cut：
  - **source 側可達集合 Rs**：`v_in ∈ Rs 且 v_out ∉ Rs` 的節點 v 即為被 cut 的節點。
  - 同一個值可能有多個 min-cut。提供選項：`--cut-side {source, sink, balanced}`
    - source：越靠近 PI 的最小 cut
    - sink：越靠近 PO 的最小 cut（對 sink 側反向可達集合求）
    - balanced：取兩者之間，使兩側邏輯面積接近（可延後實作）
- 結果：stage 0 = `v_out ∈ Rs 的節點與被 cut 的節點`，stage 1 = 其餘。

### 4.3 這裡的重點

- **反向 ∞ 邊不能漏**，否則會得到不合法（unbalanced）的 cut，而且 cut 值會偏小。
- **反向邊要接到 `u_in`，不是 `u_out`**。`v_in → u_out` 會禁止「u 被 register、同時又直接餵給 stage 0 的 fanout」
  （register 由同一節點的所有 fanout 共用，這是合法且常見的最佳解），導致 cut 值偏大。
  例：`tests/fanout_share.aag` 正確為 2、接 `u_out` 得 3；EPFL `max` 130 vs 512。
  `aig_mincut --spec-graph` 可重現舊建圖做比較。
- 不需要 `v_out → v_in` 邊：`v_out` 只能經由 `v_in` 到達，source/sink 側取出的 cut 都自然滿足一致性。
- 圖大小：2|V| 個點，約 3|E| + |V| 條邊，百萬節點可處理。

## 5. 階段 B：固定 S 的多 stage（精確解）

### 5.1 LP 形式

變數：每個節點 v 的 `s(v)`（整數），每個節點 u 的鏡像變數 `m(u)`。

```
minimize   Σ_u  W(u) · ( m(u) − s(u) )

subject to
  s(PI)  = 0                                  (所有 PI)
  s(u) ≤ s(v)                                 (每條 AIG 邊 u→v)
  m(u) ≥ s(v)                                 (每條 AIG 邊 u→v)
  m(u) ≥ S − 1                                (u 是 PO driver)
  m(u) ≥ s(u)
  0 ≤ s(v) ≤ S − 1
```

- 限制式全是差分約束 `x − y ≤ c`，矩陣為 totally unimodular，LP 最佳解是整數。
- 此 LP 的對偶是 **min-cost flow**。這就是 Leiserson–Saxe retiming 的 min-area 形式（用 mirror vertex 處理 fanout 共用）。
- 檢查：S=2 時，結果的 register 數必須等於階段 A 的 min-cut 值。

### 5.2 實作路線

1. **原型**：用 `scipy.optimize.linprog`（HiGHS）直接解 LP，驗證正確性（中小電路，約 10⁵ 變數內）。
2. **正式版**：轉成 min-cost flow，用 network simplex（C++ 可用 LEMON 或 OR-Tools `SimpleMinCostFlow`，或自己寫 cost scaling）。
3. 解出 s(v) 與 m(u) 後，每個節點 u 的 register chain 長度為 `m(u) − s(u)`。

### 5.3 邊界情況

- 若 `m(u) − s(u)` 對 PI 大於 0：代表輸入需要被延遲到後面 stage 才使用（balancing register，即 F 訊號的延遲），要在報表中另外列出數量。
- Dangling 節點（沒有 fanout 也不是 PO driver）：先移除，不進入 LP。
- Constant：不參與。

## 6. 面積權重

- `--reg-weight 10`：register 面積，以 NAND2 為單位。W(u) 預設對所有節點相同。
- `--mode count`：W = 1，只最小化 register 個數。
- 預留介面：W(u) 之後可以依位元寬度或不同 register 型別改變。
- 報表同時給出 register 個數與 `weight × 個數` 的面積。

## 7. 階段 C：延遲限制 T（可選，第二版再做）

延遲模型：每個 AND 節點單位延遲（level），之後可替換成 cell delay。

1. **Level 窗口（必要條件，直接加入 LP）**
   - `L(v)` = AIG level（PI 為 0）。
   - 每個 stage 組合深度 ≤ T，因此 `s(v) ≥ ceil(L(v)/T) − 1`（ASAP 下界）。
   - ALAP 上界由 PO 方向的 reverse level 同理推導。
   - 把這些當作 `s(v)` 的上下界加入 LP。
2. **Lazy constraint generation（充分條件）**
   - 求解後，對每個 stage 計算 stage 內最長路徑。
   - 若超過 T，找出違規路徑的起點 u 與終點 v，加入約束 `s(v) − s(u) ≥ 1`，重解。
   - 重複直到無違規，或判定不可行（報告該 T 不可行）。
3. 不要一開始就建 W/D 矩陣（O(n²)，大電路不可行）。

## 8. 輸出與資料結構

- 輸入：AIGER（`.aag` 與 `.aig`），指令列參數 `--stages S`、`--reg-weight`、`--cut-side`、`--period T`。
- 輸出：
  - `stage_assign.json`：每個節點的 `s(v)`、`m(v)`。
  - `cuts.json`：每個 cut 邊界 i（i = 1..S−1）上的訊號集合（節點 id、是否為 balancing register、chain 位置）。**這是後續 Curtis 重編碼階段的輸入。**
  - 報表：每個 cut 的寬度 R_i、總 register 數、總 register 面積、balancing register 數、最大 stage 深度。
- 內部資料結構：節點陣列 + CSR 格式的 fanin/fanout，避免用物件指標。

## 9. 驗證計畫

1. **Balanced 檢查**：DP 驗證每條 PI→PO 路徑穿過的 register 數恰為 S−1（對每個節點驗證所有 fanin 路徑的 stage 一致）。
2. **功能等價**：把 register 視為 wire（刪除後）電路必須與原 AIG 完全相同；另外可輸出插入 register 的 BLIF/Verilog，用 ABC 的 sequential equivalence 或帶 latency 的模擬比對。
3. **小電路暴力解**：節點數 ≤ 20 時，枚舉所有合法 s(v) 指派，確認 LP / min-cost flow 結果為最優。
4. **交叉驗證**：S=2 時 max-flow 結果 = LP 結果。
5. **單調性**：S 增大時，總 register 數不會減少（無延遲限制時）。
6. **測試電路**：手寫小例子（加法器、比較器、parity tree）、隨機 AIG、ISCAS / EPFL benchmark。
7. **規模測試**：EPFL 大電路（如 `hyp`、`multiplier`、`log2`）記錄執行時間與記憶體。

## 10. 里程碑

1. AIGER parser、AIG 資料結構、level 計算。
2. 階段 A（S=2 max-flow）+ 驗證 1、2、3。
3. 階段 B（先 scipy LP，再 min-cost flow）+ 驗證 4、5。
4. 輸出 `cuts.json` 與報表。
5. 階段 C（T 限制）。
6. 效能最佳化，跑大型 benchmark。

## 11. 已知限制（先寫進文件，避免誤解）

- 本階段只回答「用現有 AIG 節點當 register 位置，最少需要幾個」，這是**結構上界**，通常高於函數性下界 `⌈log₂ μ⌉ + |F_delay|`。
- 不做邏輯複製與重編碼，後續階段（SAT 刪冗餘訊號、Curtis 重編碼）才會處理。
- 延遲限制的 lazy 方法是啟發式，不保證在 T 緊的情況下找到最優解。
- 反相邊假設免費；若目標工藝庫對反相器有面積成本，要另外處理。

## 12. 待確認的設計選擇

- 實作語言（建議：Python 原型驗證 → C++ 正式版）。
- PI 是否都在 stage 0 到達，或可有不同到達 stage。
- 延遲模型使用單位延遲或 cell library delay。
- 是否需要輸出可直接給 ABC 使用的格式（BLIF / Verilog / AIGER with latches）。