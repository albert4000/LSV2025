# 計劃：Window-based random re-encoding 跳出 local minimum，推進 area/delay Pareto front

## 0. 給 Claude Code 的說明

這份計劃要回答一個問題：**把 AIG window 內部訊號所代表的 output equivalence class 做隨機重編碼，當作擾動（kick），
再交給 ABC 既有的最佳化 script 收斂，能不能在相同 runtime 下得到比 ABC 自身（含 `&deepsyn`）更好的 area/delay
Pareto front？**

工作分成五個階段。**Phase 2 與 Phase 3 結束後各有 go/no-go 判斷**，沒過門檻就停下來寫結論，不要為了讓後續階段跑起來而
調整門檻。每個階段結束時停下來回報結果再繼續。

開始前請先讀 `src/ext-lsv/lsvPipeAcd.cpp`（`findWindow`、`analyzeGroup`、`assignCode`、`synthesize`、
`decoderTruths`、`checkTruths`、`cecGia`），本計劃大部分元件從那裡改寫。下列假設任何一項不符就先回報。

### 假設（需確認）

| 項目 | 假設 |
|---|---|
| 表示法 | 組合 AIG（`Gia_Man_t`），不含 register |
| Benchmark | EPFL combinational suite，位於 `../testcases/aig/epfl/*.aig` |
| 成本模型 M1（內迴圈） | area = AND 數，delay = AIG level（unit delay） |
| 成本模型 M2（最終評估） | `read_lib ../Downloads/gradmap_libs/asap7.lib; [&dch -f;] &nf -p -R r; &put; topo; stime -p` 的 Area 與 Delay（ps），以 delay relaxation `-R` ∈ {0, 10, 25, 50, 100} 掃 trade-off。**`&nf -D` 沒有接到 mapper（設的是 `DelayTarget`，mapper 只讀 `MapDelayTarget`），不能用。** `&nf` 自己算的 Delay 與 `stime` 不一致，front 一律用 `stime` 的值 |
| Re-opt script R | 預設 `dc2`（AIG level 會變）；delay 保守版 `rs -K 8; rw; rf`（level-preserving） |
| 預算 | 每 (電路, delay 目標) 固定 wall-clock 預算，預設 小電路 60 s、大電路 600 s |
| 比較基準 | 相同 wall-clock 下的 iterated R、`&deepsyn -T`（多個 `-S` seed） |

### 名詞

- **window W = (L, C, G)**：leaves L（|L| ≤ K，由 root 往上做 reconvergence-driven 展開）、cone C（L 的 forward
  closure：所有 fanin 都在 L ∪ C 內的 AND，上限 500 個；所以 C 可以包含不通往 root 的節點）、輸出 G = C 中被 C 外部
  或 PO 使用的節點（**輸出邊界**，不是任意子集；否則舊邏輯無法刪除）。因為 G 涵蓋所有對外的使用，替換後 C 整個
  可以刪除，不需要另外算 MFFC。
- **合法性**：leaf 不可以依賴 C（否則 E 讀 leaf、leaf 又讀 G，替換後形成 loop）。C 中被任何 leaf 依賴的節點、以及 C
  內依賴它們的節點，都從 C 移除；root 被移除時 window 作廢（記為 cyclic）。這條規則是在 Phase 2 實作替換時發現的，
  Phase 1 已用此規則重跑，舊結果保留在 `results/phase1_v1_nocyclecheck/`。
- **F_W**：L → G 的局部函數；**Img(W)**：L 取遍所有值時 G 的相異向量集合；N = |Img|，r = ⌈log₂N⌉。
- **class**：F_W 值相同的 L pattern 所成的等價類（即 `analyzeGroup` 的 `inv`）。
- **code φ**：N 個 class → B^r 的 injection（固定 m = r，見第 2 節）。**encoder** E = φ∘F_W：L → B^r；
  **decoder** D：B^r → G，未使用的 code 為 don't care。正確性：D∘E = F_W 對所有 L pattern 成立，由構造保證。
- 符號：**K 一律指 window 的 leaves 上限**，class 數一律用 **N**，不要混用。
- **gain(W)** = |MFFC(C)| − AND(E) − AND(D)（只算真的會被刪除的舊 AND）。
- **move**：對一個 window 套用一個 code。**local move** 只接受 gain > 0 且 level 不違反目標；**kick** 允許 gain < 0。
- **Pareto archive**：所有已驗證（cec 通過）的 (area, delay) 非支配點，附產生它的 AIG 檔。
- **undo rate**：kick 後跑 R，結果的 (AND, level) 與 kick 前跑 R 的結果完全相同的比例。

## 1. 背景與待驗證的假說

先前在 pipeline cut 上做的確定性 ACD 重編碼（`src/develop/acd_reencode/results/phase2_results.csv`）淨面積改善只有
0–0.35%，主因是 encoder/decoder 的成本吃掉了節省量。所以本計劃**不預期重編碼本身是好的 move**；它的價值必須來自
「把電路推到 R 原本到不了的 basin，再由 R 收斂到更低的點」。以下皆為待驗證假說：

- **H1（直接 move）**：對單一 window 試 N 個隨機 code 取最好的，直接 gain > 0 的 window 比例很低（< 5%）。
- **H2（undo）**：R = `dc2` 時，affine（A）與小步數鄰域（W, s = 1）的 kick 大多被 R 還原（undo rate 高）；
  K 越大、擾動越強（U 或大 s）undo rate 越低。
- **H3（escape）**：kick + R 的 iterated local search，在相同 wall-clock 下比 iterated R 的面積更低（同 delay 目標）。
- **H4（vs SOTA）**：H3 的結果與 `&deepsyn` 相比，至少在部分電路上更好，或把兩者串接（`&deepsyn` 後再做本方法）
  能再改善。**這是整個計劃的價值判斷點。**
- **H5（delay 維度）**：改善集中在 Pareto front 的寬鬆 delay 端；緊 delay 端因 decoder 增加 level 而多數 move 被拒絕。

## 2. 重編碼的定義與編碼空間

### 定義

window 的內部訊號（輸出邊界 G）在 L 取遍所有值時，共出現 N 個相異向量，也就是 N 個 equivalence class。
重編碼就是把這 N 個 class 改用 r = ⌈log₂N⌉ 條線表示：

- 選一個 injection φ：{class 1..N} → B^r，每個 class 對應一個相異的 r-bit codeword；
- encoder 從 L 直接算出 r 條 code 線（從 truth table 重新合成，不是在 G 後面接邏輯，否則原 cone 仍在、結構沒變）；
- decoder 從 r 條 code 線還原 G，未使用的 2^r − N 個 codeword 為 don't care。

### 編碼空間有多大

- 單看一條 code 線：它是「哪些 class 取 1」的 indicator，即 class 集合的一個子集，有 **2^N** 種選擇。
  但 r 條線必須**合起來**能區分全部 N 個 class，所以不是任選 r 個子集都合法。
- 合法的 φ 總數是 (2^r)! / (2^r − N)!（把 N 個 class 放到 2^r 個 codeword 上，不可重複）。
- 其中只差在「code 線換順序」或「某條 code 線取反相」的 φ，結構上只差 inverter，算同一種。這個群有 r!·2^r 個元素，
  而且因為 N > 2^(r−1)，它在 φ 上的作用是 free 的，所以**實質不同的編碼數 = (2^r)! / ((2^r − N)! · r! · 2^r)**：

| N | r | 合法 φ 總數 | 實質不同 |
|---|---|---|---|
| 2 | 1 | 2 | 1 |
| 3, 4 | 2 | 24 | 3 |
| 5 | 3 | 6,720 | 140 |
| 6 | 3 | 20,160 | 420 |
| 7, 8 | 3 | 40,320 | 840 |
| 9 | 4 | 4.15×10⁹ | 1.08×10⁷ |
| 16 | 4 | 2.09×10¹³ | 5.45×10¹⁰ |

由此得到兩個設計限制：

- **N ≤ 4 的 window 不值得做**：只有 3 種，而且其中一種通常就是原本的結構。
- **N ≥ 9 不能窮舉，只能抽樣**；N ∈ [5, 8] 可以窮舉全部實質不同的編碼（≤ 840 個），Phase 1 用它量「隨機抽樣離最佳
  編碼多遠」。

### 抽樣方式（都在上述空間內）

| 代號 | 方式 | 用途 |
|---|---|---|
| U | 均勻隨機 φ | 最強擾動，encoder/decoder 成本最高 |
| A | affine：φ = A·c ⊕ b，c 為 `assignCode` 的 code，A 隨機可逆但**不是**排列矩陣 | 便宜的子空間 |
| W(s) | 鄰域：從 `assignCode` 的 code 出發，做 s 次「交換兩個 class 的 codeword」或「把一個 class 移到未使用的 codeword」 | 擾動強度可調；也能在 code 空間裡做 hill-climbing |
| X | N ∈ [5, 8] 時窮舉全部實質不同的編碼 | 上界與 U/A/W 的對照 |

抽樣時要先把 φ 正規化（例如固定 class 1 → 00…0、code 線依某個順序排序）再去重，避免把同構的編碼算成不同的樣本。

每個 code 的 encoder/decoder 用 `synthesize`（ISOP × don't-care fill × {raw, dc2, syn2}）合成，取 AND 最少與
level 最小的兩個候選，皆須通過 `checkTruths`。

另外加一個不靠 class 的對照擾動 **P0**：對同一個 window 直接從 truth table 重新合成 F_W（不經過 code 瓶頸）。
L 被當成自由變數，F_W 沒有 don't care，所以 P0 是確定性的（只有 1 個樣本）。若 P0 的效果與 U/A/W 相同，代表增益不是來自重編碼，而只是「局部重新合成」，結論必須如實寫出。

## 3. Phase 0：Baseline Pareto 與雜訊

1. 對每個電路產生 baseline 點（M1 與 M2 都要）：
   - iterated `dc2` 直到 AND 數不再下降；`resyn2` 同；`compress2rs` 同；
   - `&deepsyn -T t -S s`，t ∈ {60, 600}，s ∈ {1, 2, 3}；
   - delay 端：`&syn2`、`&dch -f`、`balance` 等 level 導向 script；
   - M2：上述每個 AIG 跑 `&nf -p` 與 `&dch -f; &nf -p` 兩種 flow，各以 `-R` ∈ {0, 10, 25, 50, 100} 掃 delay，用 `stime -p` 取 Area/Delay。
   - AIG 一律用 `write_aiger -s` 保留 I/O 名稱；cec 用 `cec -n`（依順序對應 I/O）。
   - Driver：`src/develop/window_reencode/phase0.py`。
2. 取所有點的非支配集合，作為 **baseline front B**。
3. 量 `&deepsyn` 不同 seed 之間的面積差異，作為雜訊尺度 σ。之後所有改善量都要和 σ 對照。

**完成條件**：`phase0_baseline.csv`（每點：circuit, script, seed, time, ands, levels, map_area, map_delay, cec）
與每電路的 front B。在此停下回報。

## 4. Phase 1：Window 統計（不改電路）

電路（依 Phase 0 結果調整）：`&deepsyn` 卡住的 B 群 adder、bar、log2、multiplier、sqrt、square、voter、sin、max、
priority，加上 A 群對照 int2float、cavlc。起點是 Phase 0 中 cec 通過、AND 數最少的 AIG（也就是已在 local minimum
的電路）。Driver：`src/develop/window_reencode/phase1.py`；指令：`lsv_reencode -a`（`src/ext-lsv/lsvReencode.cpp`）。

1. 隨機抽 1000 個 root（電路小則全取），以 reconvergence-driven 方式長出 window，K ∈ {6, 8, 10, 12}，G 取輸出邊界。
2. 統計 |L|、|C|、|MFFC|、|G|、N、r，以及 `|G| > r`（可壓縮）與 `N ≥ 5`（隨機空間夠大）的比例。
3. 對每個 window 試 U、A、W(1)、W(4)、P0 各 N_try = 16 個 code（N ∈ [5, 8] 另跑 X），記錄 gain 與 level 變化。

**完成條件**：`phase1_windows.csv`（每 window × 抽樣方式一列），附 H1 的結論。

## 5. Phase 2：單步 kick 實驗（go/no-go #1）

隨機選 window + 抽樣方式做一次 kick，然後跑 R，與「不 kick 直接跑 R」比較。每電路每種抽樣方式 200 次。

- 實作：`lsv_reencode -k`（替換 C、structural hashing、每次 kick 都對 x 做 cec），driver `phase2.py`。
- 每次 kick 都從同一個 x 出發；只挑 N ≥ 5 的 window；E 與 D 用 min-AND 候選。
- R = `&put; dc2; &get -n`（與 Phase 0 的 `dc2` 相同）。K ∈ {6, 8}（Phase 1 顯示 K ≥ 10 時多數隨機編碼超出預算）。
- 對照組：base（確定性 ACD 編碼）與 P0（不經過 code 瓶頸的直接重新合成）。若對照組的效果與 U/A/W 相同，代表
  增益來自「局部重新合成」而非重編碼。
- noop：kick 本身沒有改變 (AND, level)；undo 同時統計全部 kick 與排除 noop 後的比例。

量測：undo rate、Δarea 分布（R 之後）、Δlevel 分布、最好的 Δarea。

**Go 條件**（全部成立才進 Phase 3）：

- 至少一種抽樣方式的 undo rate < 80%；
- 該抽樣方式在 ≥ 3 個電路上出現 R 之後 Δarea < 0 的情況（比原本的 local minimum 更好），且最佳改善 > max(σ, 0.5%)。
  （R 是確定性的，所以 Δarea < 0 本身就是真實的改善；門檻是為了排除只省 1、2 個 AND 的情況。兩種判定都要回報。）

**完成條件**：`phase2_kicks.csv`，附 H2 的結論。在此停下回報。

### Phase 2 結果（`results/phase2/`，28,664 次 kick，cec 28,663 ok、1 undecided）

- **判定：NO-GO（嚴格）**。所有抽樣方式的 undo 都 < 80%（重編碼 12–18%，對照組 25–34%），且都在 3–5 個電路上
  有 Δarea < 0（寬鬆判定成立）；但沒有任何電路的最佳改善超過 max(σ, 0.5%)。最佳 Δarea：log2 −6（門檻 143）、
  sin −3（24）、sqrt −2（91）、multiplier −2（122）、priority −2（2.0）、int2float −1、cavlc −1，其餘 ≥ 0。
- **H2 不成立，方向相反**：kick 不是被 `dc2` 還原，而是 `dc2` 處理不掉——重編碼 kick 有 73–81% 在 R 之後仍大於
  R(x)（kick 本身 AND 中位數 +22）。例外是 bar（undo ≈ 0.9）。
- **增益不是重編碼特有的**：P0（直接重新合成）的變好比例 3.4–4.1% 高於 U 的 1.8–2.5%，最佳值相同。
- 實作期間發現並修正：cyclic window（見名詞）；ABC Fraig CEC 在部分 sqrt miter 上 assertion 失敗
  （`Dar_ManRwsat → Aig_ManDupDfs`），改在 fork 的子程序中跑，失敗時退回 SAT sweeping（`cecGiaSafe`）。

### Phase 2b（診斷，計劃外）：更強的 R

目的：區分「R 太弱」與「重編碼沒有額外價值」。在 Phase 2 有逃脫的 log2、sin、sqrt、priority 上，K = 8、六種
抽樣方式各 200 次，R 換成 `dc2` ×5 與 `compress2rs` ×3（含 `resub -K 6..12 -N 2`）。一輪 `&deepsyn`
（`&dch -f; &if -a -K 6; &mfs; compress2rs`）在起點上比 x 差 4–6%、單次最長 65 s，屬於重啟而不是局部最佳化，
不採用。Driver：`phase2b.py`，結果在 `results/phase2b/<R>/`。

結果（各 4800 次 kick，cec 4798 ok、2 undecided）：

- `compress2rs` ×3 明顯較能吸收 kick：重編碼 undo 38–48%（`dc2` 時 12–18%）；log2、sin 上 R 之後 Δarea 中位數
  為 0，變好比例 18–50%。`dc2` ×5 則幾乎吸收不了（log2 變差 98–100%）。
- 但改善量沒有變大：最佳 Δarea log2 −12、sin −7、sqrt −2、priority −2，全部低於門檻，嚴格判定仍是 NO-GO。
- 和 R 的不動點比較：對 x 反覆跑 `compress2rs` ×3 會收斂到 log2 28495、sin 4866，所以單次 kick 真正超越
  不動點的只有 log2 −8、sin −7（≤ 0.15%）。`dc2` ×5 反覆執行 50 次仍持續下降（log2 −61），它的「逃脫」與
  未收斂無法區分。
- 對照組 base/P0 的最佳值與變好比例都 ≥ 重編碼（log2 P0 變好 49.5%，U 25%）。
- **結論**：Phase 2 的瓶頸確實部分是 R 太弱（吸收率），但即使 R 夠強，單次 kick 的增益仍在 0.1% 以下，而且
  不是重編碼特有的。

### Phase 2c（診斷，計劃外）：未最佳化的起點與結構變化

問題：重編碼能不能大幅改變電路結構？改變的結構有什麼特徵？

- 起點 x0 = 原始 EPFL AIG（不先 `&deepsyn`），18 個電路（Phase 1 的 12 個加 arbiter、mem_ctrl、i2c、ctrl、
  router、dec）。R = `compress2rs` ×3。抽樣 U、A、W4、base、P0，K = 8。
- kick 規模：m1（1 個 window，200 次）；P10 / P50（在前一個結果上連續挑 window 重編碼，直到被替換的 cone 合計
  達 AND(x0) 的 10% / 50%，各 20 次）。實作：`lsv_reencode -k -m / -P / -W`。
- 結構量測（`structdiff.py`）：2048 個隨機 pattern 模擬，節點函數（不分正反相）作為 key。new_frac = R(kick) 中
  函數不在 R(x0) 出現的 AND 比例；kept_frac = R(x0) 的 AND 在 R(kick) 仍存在的比例；另記 XOR/MUX 數、共用率、
  新節點的平均相對 level。參考尺度：x0、R*(x0)、Phase 0 最佳 AIG 對 R(x0) 的同一組量測。
- Driver：`phase2c.py`，結果在 `results/phase2c/`。

結果（21,446 次 kick，cec 全部通過，1 次 undecided）：

- **結構可以大幅改變，但改變後面積變差，R 也救不回來**。P50 的 new_frac 中位數 U 0.44、W4 0.47，與 `&deepsyn`
  相對 R(x0) 的改變量（0.13–0.86）同一個量級；但 R 之後 AND 數中位數 U +48%、W4 +59%、A +37%、base +21%，
  P0 只有 +1%。P10 時 U +12%、W4 +14%。同樣幅度的結構改變，`&deepsyn` 會變好、重編碼會變差，所以問題不在改變的
  多寡，而在方向。
- 懲罰大小的排序在所有電路上都一致：W4 > U > A > base > P0，也就是 code 越隨機懲罰越大。A（affine）會讓 XOR 數
  在 P50 增加為約 3.9 倍（code bit 是 base bit 的 XOR）。
- **P0 連續套用本身就是有效的最佳化**：mem_ctrl P50 −8.8%、P10 −6.9%，bar −3.6～−4.3%，arbiter −0.75%，
  而且 100% 的 kick 都變好。這是「8-input 多輸出 window 的 exact 重新合成」，與重編碼無關。
- 單一 window（m1）：結構改變很小（new_frac 中位數 0.2–2%）。18 個電路中有 16 個的 new_frac 與 d_R 是
  正相關（改變越多、結果越差；Spearman 0.11–0.99）；例外是 bar（−0.59）與 arbiter（−0.08）。
- 與 R 的不動點 R*（`compress2rs` 反覆 12 次）比較後，真正超越 R* 的單次 kick 增益：bar −177（−5.6%）、priority −14、
  router −14、ctrl −11、int2float −8、cavlc −8、i2c −6；log2、sin、mem_ctrl 的「改善」都還不到 R*，算術電路
  （multiplier、square、sqrt、voter、adder）只有 0 到 −1。`&deepsyn` 在這些電路上是 −13% 到 −84%。
- 只有 router 上重編碼有對照組沒有的增益：U −9.7%、A −6.2%、W4 −4.1%，而 base −2.1%、P0 0。最佳改善同時
  降低 level（router 20 → 16，priority 60 → 34）；觸發它們的都是小 window（|L| = 3–8，N = 5–18），kick 本身
  只多 6–64 個 AND。
- 很多不同的 root 會收斂到同一個 AND 數，例如 arbiter 有 154 個 root 都得到 11838，bar 有 62 個 root 都得到 3137。
  也就是說，逃脫主要決定於 R 的少數幾個吸引點，與 kick 的位置關係不大。

### Phase 2d（診斷，計劃外）：放大 window 加上隨機編號

`phase2d.py`：Phase 2c 的 m1 設定，K ∈ {6, 8, 10, 12}，抽樣方式 U 與 P0，8 個有逃脫的電路，各 200 次 kick。

- **放大 K 並不會讓 U 真的用到大 window**：K = 12 時，U 的 kick 只有 20% 用到 |L| = 12，48% 的 |L| 只有 3–4
  （P0 則有 41% 用到 |L| = 12）。大 window 的隨機編碼大多超出預算（Phase 1 中 U 的 fail rate 在 K = 10、12 時
  為 67–69%），能通過的大 window 也都是 N 小的（中位數 9）。
- 能通過的大 window 結果比較差：|L| 9–12 的 kick 有 14% 變好、73% 變差；|L| 3–4 則是 41% 變好、44% 變差。
- 各電路的最佳值幾乎不隨 K 改變（router U 在 K = 6–10 都是 −14，K = 12 時是 −5）；P0 在 K ≥ 10 時追上
  （router −13）。

**放寬預算**（`lsv_reencode -B b`：每個合成電路的上限為 b·(|C| + 10)，預設 b = 3，0 表示不設上限；只跑 U，
`phase2d.py --budget 10` 與 `--budget 0 --nocec`。B = 0 時 kick 後的 AIG 會膨脹數倍，cec 要數百秒，所以跳過每次
kick 的全電路 cec；每個 window 仍然窮舉驗證 E∘D = G，每個 job 的最佳 R(kick) 也都通過 cec（32/32）。）

K ∈ {10, 12} 合併，依實際 |L| 分組：

| 預算 | \|L\| ≥ 9 的比例 | \|L\| 9–12 的 N 中位數 | kick 後 AND 增加（中位數） | R 之後變好 | R 之後 ΔAND 中位數 | 最佳 |
|---|---|---|---|---|---|---|
| 3(\|C\|+10) | 25% | 9 | +1.3% | 14% | +0.1% | −4.9% |
| 10(\|C\|+10) | 30% | 18 | +4.8% | 6% | +4.4% | −1.6% |
| 不設上限 | 48% | 84 | +53% | 2% | **+67%** | −1.6% |

- 放寬預算後，大 window 的隨機編碼確實能通過，但電路大幅膨脹（router K = 12 時 kick 後中位數 +11,933 AND，
  原電路只有 257），`compress2rs` 幾乎完全救不回來。
- 各預算下的最佳結果都來自 |L| = 3–4 的 window（−12.2%，ctrl），與預算無關。放寬預算只會讓最佳值持平或變差：
  router K = 8 從 −9.7% 變成 −5.5%，K = 12 從 −3.4% 變成 0；priority K = 12 從 −3.2% 變成 −0.9%。
- 有些電路即使不設上限，|L| 也到不了 K（ctrl 在 K ≥ 8 時為 0%）：這是電路本身的限制（findWindow 碰到 PI，或
  cone 超過上限），不是預算造成的。

### `&deepsyn` 元件 ablation（原始 EPFL、各步驟後接 `compress2rs`、迭代 8 次、K = 6）

| 電路 | rs | lut | dch | mfs | full | Phase 0 最佳 |
|---|---|---|---|---|---|---|
| mem_ctrl | 38944 | 29645 | 26040 | 14090 | 12472 | 7871 |
| router | 145 | 164 | 138 | 108 | 108 | 111 |
| i2c | 1016 | 1003 | 984 | 966 | 943 | 684 |
| int2float | 207 | 205 | 189 | 190 | 177 | 121 |
| priority | 428 | 433 | 406 | 422 | 400 | 374 |
| bar | 3141 | 2688 | 2688 | 2688 | 2688 | 2688 |
| arbiter | 11839 | 11838 | 11839 | 11838 | 11839 | 1920 |
| sin | 4972 | 5233 | 5087 | 5194 | 5186 | 4888 |

lut = `&if -a -K 6; &st`；dch 在前面加 `&dch -f`；mfs 是在 `&if` 之後加 `&mfs -e -W 20 -L 20`；full 兩者都加。arbiter
只有加上 `&fx` 才會動（`&if; &mfs; &fx` 三輪後為 3507），或是做 `collapse`（2198）。

結論：`&deepsyn` 的主要動力是 `&mfs`（以 SAT 計算 window 內的 SDC/ODC 做 don't-care 最佳化），其次是 LUT 化之後
重新拆解（bar）、choices（`&dch`），以及跨節點的代數抽取（`&fx`，arbiter）。另外，它每一輪都從目前的結果繼續，
即使變差也不退回（random walk），只記錄最佳解。這些機制在重編碼加 `compress2rs` 裡都沒有：window 沒有用到
don't care、E/D 不和 window 外共用、每次 kick 都從同一個 x 重新開始。

### Phase 2e（診斷，計劃外）：null kick 與殘留分析

**指標**：所有 Phase 2 系列的 d_R = AND(R(kick(x0))) − AND(R(x0))，也就是 AIG 的 2-input AND 數（不是 mapping 後的
面積）；baseline 是同一個 R 套在未 kick 的起點 x0 上。

**null kick（抽樣 N，`lsv_reencode -M N`）**：把 window 原封不動放回去，AIG 完全相同，只有節點編號順序和真正的
kick 一樣被打亂。R 的啟發式依賴節點順序，所以 N 量的是純順序雜訊。`phase2e.py`：Phase 2c m1 的設定（原始起點、
`compress2rs` ×3、K = 8、200 次，seed 相同，所以 U … P0 的 kick 與 Phase 2c 相同），加上 N；Phase 2 的設定
（`&deepsyn` 起點、`dc2`）也補跑 N（`results/phase2/raw/*_N_K*.csv`）。

- 順序雜訊因電路而異：adder、bar、multiplier、sqrt、arbiter、priority、ctrl、router、dec 完全沒有；
  log2 −11～+13、mem_ctrl −94～+91、i2c −4～+13、sin −4～+8、cavlc −3～+5、max −2～0、voter −1～0。
- **修正 Phase 2 / 2b 的結論**：Phase 2 設定下，log2 的 N 範圍是 −7～+9、sin 是 −3～+2，而各抽樣方式在 log2 的最佳值
  是 −3～−6、在 sin 是 −3，全部落在雜訊範圍內。Phase 2 中 log2、sin 的「逃脫」只是節點順序造成的雜訊。
- 原始起點下，超出雜訊的最佳改善（重編碼最佳 / 對照組 P0、base 最佳 / N 最佳）：bar −181 / −176 / 0、
  mem_ctrl −719 / −730 / −94、priority −18 / −9 / 0、router −14 / −3 / 0、i2c −21 / −19 / −4、ctrl −11 / −10 / 0、
  cavlc −13 / −17 / −3、sin −11 / −12 / −4、int2float −8 / −7 / 0；log2 −17 / −13 / −11 與 max −2 / −2 / −2
  沒有超出雜訊。重編碼明顯優於對照組的只有 router 與 priority。

**殘留分析（`residue.py`；`lsv_reencode -W` 記錄 kick 後每個節點的來源：X 原本、E、D、T 重建的 TFO）**：以函數
（模擬簽章，不分正反相）比對 R(y) 與 R(x)。在沒有順序雜訊的電路上，U/A/W4 的 d_R 可拆成：

| | d_R | E/D 殘留 | x 的函數（保留 − 失去） | 新函數（window 下游） | 新函數（其他位置） |
|---|---|---|---|---|---|
| U 變差 | +43.3 | **+34.0** | −10.2 | +16.3 | +1.6 |
| U 變好 | −23.8 | +1.4 | **−84.6** | +28.1 | +31.7 |

- **變差的主因是 E/D 沒被消掉**（約 80%）。關鍵指標是「新的 code 線」（不是 x 原本就有的函數）是否還存在於 R(y)：
  只要有一條留下來，92% 變差（d_R 中位數 +28）；全部消失時 30% 變好、24% 變差（中位數 0）。
  這個關係在每個電路內都成立（有新 code 線留下時，17 個電路中有 16 個變差率 ≥ 85%）。
- E/D 的存活率隨 N 與 r 上升：N = 5 時 12%，N = 6–8 時 26%，N = 17–32 時 44%，N > 32 時 70%。r = 3 時變好 27%、
  r = 6–8 時 0%。在單一電路內，N ≤ 8 或 kick 成本 ≤ 10 也能提高變好的比例，但多數電路的最佳 kick 落在 N > 8、
  成本 > 10（priority −18、mem_ctrl −719、i2c −21）。若用 N ≤ 8 且成本 ≤ 10 過濾，無雜訊電路的變好比例從 20% 升到 37%，
  但 7 個電路中只有 2 個保留了最佳 kick。
- **XOR 型 code 線不是原因（更正）**：合併所有電路時，code 線平均 ANF degree ≤ 1.5 的 kick 存活 63%、變差 79%，
  但這類 kick 幾乎只出現在 ctrl、dec、cavlc、int2float，是電路造成的混淆。在同一電路內比較時，degree 高低對結果沒有
  一致的影響（例如 ctrl 的低 degree 組反而變好 38%、高 degree 組 31%）。抽樣方式 A（code 為 base code 的 XOR）
  與 U 相比，變差率較低的電路有 7 個、較高的 8 個、相同的 3 個。
- 變好的 kick：E/D 幾乎完全被消掉，增益來自 R(y) 不再需要許多 x 原本的函數（其中 OTHER 區佔 24.5%、TFO 佔 11.5%），
  再由 R 產生新的函數取代，而且有一半落在 window 下游之外：改變不只發生在 window 內。

### Phase 3a：ILS 對 `&deepsyn`（只看面積）

設定（`phase3a.py`、`lsv_reencode -i`）：
- 起點是原始 EPFL AIG，R = `compress2rs` ×3，K = 8，每次 kick 一個 window。
- 若 AND 變少，或 AND 相同但 level 變低，就接受；`-e` 變體也接受 AND 相同且 level 不變差的結果（平移）。
- 每個 job 20 分鐘（wall clock），seed 1；對照組為 `&deepsyn -T 1200 -S 1`。
- 每次 kick 都 cec；最終結果對 x0 做 cec，144 個 run 全部通過；超過 Phase 0 最佳值的結果另外用 `&cec` 確認。

最終 AND 數（ILS 欄取 8 種設定中最好的一個）：

| 電路 | R(x0) | ILS 最佳（設定） | P0 最佳 | N 最佳 | `&deepsyn` 20 分鐘 | Phase 0 最佳 |
|---|---|---|---|---|---|---|
| router | 145 | **97**（U-e，lev 13） | 145 | 145 | 111（lev 17） | 111 |
| priority | 432 | **367**（A，lev 96） | 399 | 427 | 376（lev 49） | 374 |
| ctrl | 90 | **66**（A） | 74 | 90 | 68 | 66 |
| bar | 3141 | **2911**（A） | 2988 | 3141 | 2952 | 2688（syn2） |
| voter | 7936 | **7825**（P0） | 7825 | 7936 | 7935 | 7932 |
| sqrt | 18253 | **18225**（P0） | 18225 | 18252 | 18263 | 18232 |
| sin | 4991 | **4879**（P0-e） | 4879 | 4957 | 4888 | 4888 |
| multiplier | 24315 | **24294**（P0-e） | 24294 | 24315 | 24321 | 24315 |
| square | 15800 | 15792（A） | 15793 | 15796 | 15798 | 15789 |
| log2 | 29077 | 28952（P0-e） | 28952 | 28978 | **28583** | 28583 |
| max | 2824 | 2779（U-e） | 2822 | 2823 | **2677** | 2677 |
| i2c | 1030 | 828（A-e） | 883 | 1015 | **684** | 684 |
| int2float | 207 | 144（A-e） | 178 | 207 | **122** | 121 |
| cavlc | 623 | 526（U-e） | 546 | 621 | **375** | 382 |
| arbiter | 11839 | 8790（U-e） | 11722 | 11839 | **1998** | 1920 |
| mem_ctrl | 41030 | 35776（U） | 36392 | 38377 | **7582** | 7871 |
| adder / dec | 892 / 304 | 不變 | | | 不變 | |

發現：
- **ILS 累積了單次 kick 的增益**：router 單次 kick 最佳 −14，ILS 達到 −48；priority −18 → −65；ctrl −11 → −24；bar −181 → −230。
- **9 個電路勝過 `&deepsyn`（20 分鐘）**，6 個超過 Phase 0 的所有變體：router（−12.6%）、priority、voter（−1.3%）、sqrt、sin、multiplier。
  - 控制類電路（router、priority、ctrl、bar）靠重編碼（U/A），P0 遠遠落後。
  - 算術類電路（voter、sqrt、sin、multiplier、log2）則是 P0 最好，code 瓶頸反而沒有幫助。
- **只打亂節點順序（N）幾乎沒有進展**：無雜訊電路完全不動；mem_ctrl 的 −2653 和 log2 的 −99 是順序雜訊累積的結果。所以上述增益確實來自 kick。
- **平移（-e）利弊參半**：
  - 有幫助：router U 從 114 降到 97、arbiter 從 9766 降到 8790、int2float A 從 173 降到 144。
  - 有害：mem_ctrl、priority U。
- **`&deepsyn` 大勝的電路**（arbiter、mem_ctrl、cavlc、i2c、int2float、max、log2）正是 ablation 中靠 `&mfs` 的 don't care 或 `&fx` 的電路。window 重編碼保留 G 的函數，沒有用到 don't care，所以追不上。
- **大電路受限於迭代次數**：每次 iteration 都對整個 AIG 跑 R。mem_ctrl、log2、multiplier 在 20 分鐘內只有 108～187 次 iteration，而且 arbiter、mem_ctrl 的曲線到 1200 秒仍在下降。
- **只看面積會犧牲 delay**：priority 的 level 從 `&deepsyn` 的 49 變成 96。

### Phase 3a 多 seed 與 Phase 3b：ILS 疊在 `&deepsyn` 之上

> **注意（更正）**：本節與 Phase 3a 表格中的「ILS 最佳」是 6～8 種 ILS 設定（各 1200 s）取最小值，而 `&deepsyn` 只跑
> 一次 1200 s，因此 ILS 等於用了數倍的運算量。改用單一固定設定比較時，A-e 在原始起點下有 8 個電路贏過 `&deepsyn`，
> A / A-e 在 `&deepsyn` 之後有 11 個電路贏過繼續跑的 `&deepsyn`。詳見 `window_reencode_experiments.md` 第 4 節。

設定：
- Phase 3a 補跑 seed 2、3。設定為 U、A、P0 及其 -e 變體，加上 `&deepsyn`；排除 adder 和 dec。
- Phase 3b 的起點是 Phase 3a 中 `&deepsyn` seed 1 的結果。
  - ILS 用 seed 2，設定為 U、A、P0、N 及 U/A/P0 的 -e 變體。
  - 對照組是從同一起點繼續跑 `&deepsyn -S 2`。
  - 每個 job 都是 20 分鐘。
- `lsv_reencode -i` 的起點改為 x 與 R(x) 中較好的一個；原始起點下兩者等於 R(x)，所以 Phase 3a 不受影響。
- cec 結果：Phase 3b 的 112 個 ILS run、seed 2/3 各 96 個 run 全部通過；關鍵結果另外用 `&cec` 確認。

Phase 3a 的 3 個 seed（ILS 欄是每個 seed 下 6 種設定的最小值）：

| 電路 | ILS S1/S2/S3 | `&deepsyn` S1/S2/S3 | ILS 勝場 |
|---|---|---|---|
| router | 97/97/96 | 111/112/111 | 3/3 |
| priority | 367/367/366 | 376/378/374 | 3/3 |
| bar | 2911/2923/2920 | 2952/2952/2952 | 3/3 |
| voter | 7825/7844/7849 | 7935/7934/7932 | 3/3 |
| sqrt | 18225/18224/18229 | 18263/18232/18265 | 3/3 |
| sin | 4879/4902/4884 | 4888/4929/4920 | 3/3 |
| multiplier | 24294/24298/24295 | 24321/24367/24334 | 3/3 |
| square | 15792/15790/15791 | 15798/15819/15789 | 2/3 |
| ctrl | 66/69/68 | 68/67/66 | 1/3 |
| log2、max、i2c、int2float、cavlc、arbiter、mem_ctrl | | | 0/3 |

- seed 1 的結論穩定：7 個電路 3 個 seed 全勝，ctrl 則勝負互見。
- 若所有電路固定用同一種設定，A-e 最好：48 場中贏 21 場；但以幾何平均計，AND 數仍是 `&deepsyn` 的 1.28 倍，主要被 arbiter 和 mem_ctrl 拉高。

Phase 3b（起點是 `&deepsyn` seed 1 的結果）：

| 電路 | 起點 | ILS 最佳（設定） | N | `&deepsyn` 再跑 20 分鐘 |
|---|---|---|---|---|
| router | 111 | **99**（U-e，−10.8%） | 111 | 111 |
| ctrl | 68 | **63**（U-e/A-e，−7.4%） | 68 | 67 |
| int2float | 122 | **117**（A-e，−4.1%） | 122 | 122 |
| cavlc | 375 | **368**（U-e，−1.9%） | 375 | 375 |
| sin | 4888 | **4798**（U，−1.8%） | 4854 | 4888 |
| priority | 376 | **370**（U-e，−1.6%） | 376 | 376 |
| bar | 2952 | **2911**（A，−1.4%） | 2952 | 2952 |
| log2 | 28583 | **28433**（U，−0.5%） | 28466 | 28583 |
| voter / sqrt / multiplier / square | | −12 / −10 / −23 / −10 | −4 / 0 / −19 / −6 | −8 / 0 / 0 / 0 |
| max | 2677 | 2676 | 2677 | 2676 |
| i2c | 684 | 678 | 684 | **675** |
| arbiter | 1998 | 1975 | 1998 | **1887** |
| mem_ctrl | 7582 | 7548 | 7575 | **7370** |

- **ILS 能在 `&deepsyn` 的結果上繼續縮小**：16 個電路全部變小（max 只少 1），而同樣多跑 20 分鐘的 `&deepsyn` 只有 6 個變小。
  - 12 個電路 ILS 勝出：ctrl 63、int2float 117、cavlc 368、sin 4798、log2 28433 都低於 Phase 0 所有變體的最佳值。
  - max 打平；i2c、arbiter、mem_ctrl 仍是繼續跑 `&deepsyn` 比較好。
- **扣除順序雜訊**：在有雜訊的電路上，N 本身也有進展（log2 −117、multiplier −19、sin −34），所以 log2、multiplier 的增益大多在雜訊範圍內。sin 的 −90 和控制類電路的增益則不是雜訊（這些電路的 N 為 0）。
- **起點的影響**：router 從原始電路出發的 ILS（96～97）比從 `&deepsyn` 出發（99）更好；ctrl、int2float、cavlc 則是從 `&deepsyn` 出發較好。
- **level**：
  - 持平或更低：ctrl（8/8）、cavlc（38→37）、log2（350→347）、priority（49→38）。
  - 變高：router（17→18）、int2float（13→16）、sin（181→187）。

### Kick 的位置與改善機率（`position.py`，Phase 2e 的單次 kick）

方法：
- 對每次 kick 的 root 計算位置與結構特徵：
  - 相對位置 rel_pos = lev / (lev + down)，0 代表靠近 PI，1 代表靠近 PO。
  - fanout、support、能到達的 PO 數。
  - 2 層內是否有 XOR / MUX 結構、2 層內的最大 fanout（近似 control signal）。
  - 以及 window 的大小。
- 在每個電路內計算 AUC，也就是「變好的 kick 的特徵值大於其他 kick 的機率」，0.5 代表沒有關係。
- null kick（N）用同樣的方式分析，作為雜訊對照。

結果：
- **靠近 PI 的 window 比較容易變好**，這在沒有順序雜訊的電路上成立：

  | 電路 | 靠近 PI | 中間 | 靠近 PO |
  |---|---|---|---|
  | bar | 93% | 14% | 8% |
  | arbiter | 29% | 0% | 0% |
  | ctrl | 65% | 41% | 27% |
  | multiplier | 5% | 0% | 0% |

  - sqrt、priority、router 沒有位置差異。
  - 控制 kick 成本之後效應仍在（bar 成本 > 10 的 kick：98% 對 10%）。所以不只是因為「靠近 PI 的 window 比較小」。
  - 解讀：靠近 PI 的 window，它的 TFO 涵蓋很多 PO；E/D 被消掉之後，R 會重新表達一大片下游電路。這和 Phase 2e「增益有一半落在 window 下游以外」的發現一致。
- **有雜訊的電路不能下這個結論**：log2、max 的 null kick 也呈現同樣的位置效應（max 靠近 PI 的 N 有 53% 變好）。靠近 PI 的 kick 會重新編號更大的 TFO，所以順序雜訊也比較大。
- **XOR / MUX 鄰近與否沒有關係**（AUC 約 0.5）；例外是 multiplier 的 XOR（0.69）。
- **鄰近的高 fanout 訊號**在 bar（0.89，barrel shifter 的 shift 控制線）、square、arbiter、sqrt 上與改善正相關，在 ctrl、priority、max 上負相關，所以沒有一致的規律。
- **window 大小對兩種方法的影響相反**：
  - 重編碼偏好小 window：無雜訊電路的 AUC，nC 0.41、N 0.36、成本 0.40。
  - P0 偏好大 cone：nC 0.83、nG 0.85。
- **電路類型**：「算術或控制」不是決定因素。bar（算術）的改善率最高，達 70%；adder、dec 都是 0%。真正的關鍵是 R 對該電路是否已經卡住，以及 window 是否位於上游。

### 局部 R（`lsv_reencode -i -L h`，初步）

設計：
- 每次 iteration 只對 kick 替換部分（E/D 節點）周圍 h 步（沿 fanin 與 fanout）的區域跑 R。
- 區域必須凸：凡是經由外部節點依賴區域的節點，都移出區域，確保優化後能整塊放回而不形成 loop。所以區域不會超過 h 步的範圍。
- 被接受的結果都會和 cur 做 cec。

初步測試（4 個並行，每個 120 秒，U-e，`-c`；括號內為 iteration 數）：

| 電路 | 全電路 R | L=2 | L=4 | L=8 |
|---|---|---|---|---|
| arbiter | 11827（36） | 9500 | **8348**（1393） | 10895 |
| bar | 2978（184） | 2932 | **2905**（639） | 2947 |
| sin | 4953（45） | 4964 | 4975 | **4948**（257） |
| router | **97**（2753） | 115 | 108 | 100 |

- arbiter 的 L=4 只跑 120 秒，就比 Phase 3a 全電路 R 跑 20 分鐘（8790）更好。
- router 這類小電路沒有加速效果，而且局部 R 失去了全域整理的能力。
- level 可能大幅升高（arbiter L=8：88 → 294）。
- 正式比較（Phase 3d，20 分鐘）結果見下方「Phase 3d」。

### Phase 3c：交替流程 `&deepsyn` → ILS（`phase3c.py`，`results/phase3c/`）

每個 job 總預算 3600 s，16 個電路 × seed 1、2，128 個 job 全部通過 `&cec`（設計見 `window_reencode_experiments.md` 第 3 節）。

- **D**：`&deepsyn -T 3600`，一次跑滿。
- **DR**：`&deepsyn -T 600` × 6，每段換 seed 從上一段結果重開。
- **DI**：(`&deepsyn -T 600`; ILS U-e 600 s) × 3。
- **DN**：同 DI，但 ILS 改用 null kick（N-e）。

最終 AND 數（粗體為該 seed 中最好的）：

| 電路 | D S1 | DR S1 | DI S1 | DN S1 | D S2 | DR S2 | DI S2 | DN S2 |
|---|---|---|---|---|---|---|---|---|
| bar | 2952 | 2952 | **2946** | 2952 | 2952 | 2952 | 2952 | 2952 |
| log2 | 28694 | 28539 | 28477 | **28473** | 29274 | 28832 | **28581** | 28630 |
| multiplier | 24352 | **24286** | 24295 | 24301 | 24985 | **24676** | 24882 | 24708 |
| sqrt | 18228 | 18225 | **18181** | 18184 | 18196 | 18196 | **18177** | 18180 |
| square | 15811 | 15807 | **15786** | 15794 | 15959 | 15909 | 15861 | **15851** |
| voter | 7969 | 7934 | **7928** | 7931 | 8327 | 8098 | **8025** | 8071 |
| sin | 4916 | 4854 | **4788** | 4850 | 4976 | 4887 | **4818** | 4907 |
| max | 2676 | 2676 | 2676 | 2680 | 2676 | 2676 | 2676 | 2676 |
| priority | 379 | 379 | **368** | 379 | 373 | 373 | **371** | 373 |
| int2float | 124 | **109** | 118 | 124 | 124 | 125 | 127 | **121** |
| cavlc | **359** | 369 | 361 | 372 | **360** | 369 | 365 | 361 |
| arbiter | 1743 | 1626 | **785**（lev 167） | 1893 | 1623 | 1549 | **1435** | 1684 |
| mem_ctrl | 7078 | **7039** | 7260 | 7352 | **7070** | 7118 | 7348 | 7435 |
| i2c | 691 | 689 | **679** | 690 | 690 | 686 | **677** | 690 |
| ctrl | 67 | 67 | **65** | 67 | 66 | 66 | **64** | 66 |
| router | 110 | 110 | **99**（lev 12） | 110 | 116 | 116 | **95**（lev 12） | 116 |

32 組（電路 × seed）兩兩比較：

| 比較 | 贏 | 平 | 輸 |
|---|---|---|---|
| DI 對 D | 24 | 3 | 5 |
| DI 對 DR | 23 | 3 | 6 |
| DI 對 DN | 25 | 2 | 5 |
| DN 對 DR | 11 | 9 | 12 |
| DR 對 D | 17 | 11 | 4 |

- 交替有效：DI 在 32 組中有 22 組是四種設定裡最好的。
- 重開本身有幫助（DR 對 D 17:4），但 DI 仍然贏 DR 23:6。
- 增益來自重新編碼的 kick：DN 和 DR 打平（11:12），DI 對 DN 25:5。
- DI 輸的電路：mem_ctrl（兩個 seed）、cavlc、int2float，和 3a 一致，這些電路上 `&deepsyn` 的時間比較值錢。
- level：`&deepsyn` 只看 AND。arbiter S1 的 785 是第 5 段 `&deepsyn` 在 ILS 結果上找到的，level 從 32 升到 167，不是 Pareto 改善。

### Phase 3c 追加：交替流程中的 ILS 改用局部 R（DLU、DLA）

和 DI 相同的分段（`&deepsyn` 600 s 與 ILS 600 s 交替 3 輪），ILS 改為局部 R：
- **DLU**：`-M U -e -c -L 4`。和 DI 只差在局部 R（以及 `-c`）。
- **DLA**：`-M A -e -c -L 4`，即 Phase 3d 最穩定的設定。

第 1 段 `&deepsyn` 和其他設定完全相同。64 個 job 全部 cec ok。最終 AND 數（`*` 為該列最佳）：

| 電路 | D S1 | DR S1 | DI S1 | DLU S1 | DLA S1 | D S2 | DR S2 | DI S2 | DLU S2 | DLA S2 |
|---|---|---|---|---|---|---|---|---|---|---|
| bar | 2952 | 2952 | 2946 | 2773 | **2720** | 2952 | 2952 | 2952 | 2790 | **2745** |
| log2 | 28694 | 28539 | 28477 | **28385** | 28455 | 29274 | 28832 | 28581 | 28608 | **28572** |
| multiplier | 24352 | 24286 | 24295 | 24266 | **24241** | 24985 | **24676** | 24882 | 24753 | 24739 |
| sqrt | 18228 | 18225 | **18181** | 18183 | 18182 | 18196 | 18196 | **18177** | 18181 | 18181 |
| square | 15811 | 15807 | 15786 | **15784** | 15786 | 15959 | 15909 | 15861 | 15831 | **15816** |
| voter | 7969 | 7934 | 7928 | 7900 | **7858** | 8327 | 8098 | 8025 | 7824 | **7750** |
| sin | 4916 | 4854 | **4788** | 4820 | 4824 | 4976 | 4887 | **4818** | 4841 | 4839 |
| max | 2676 | 2676 | 2676 | 2676 | 2676 | 2676 | 2676 | 2676 | 2676 | 2676 |
| priority | 379 | 379 | 368 | **365** | **365** | 373 | 373 | 371 | **367** | 372 |
| int2float | 124 | **109** | 118 | 119 | 120 | 124 | 125 | 127 | 117 | **115** |
| cavlc | **359** | 369 | 361 | 361 | 362 | **360** | 369 | 365 | 361 | 362 |
| arbiter | 1743 | 1626 | **785**（lev 167） | 1535（lev 38） | 1710 | 1623 | 1549 | **1435** | 1455 | 1522 |
| mem_ctrl | 7078 | **7039** | 7260 | 7251 | 7165 | **7070** | 7118 | 7348 | 7353 | 7379 |
| i2c | 691 | 689 | 679 | **672** | 675 | 690 | 686 | 677 | 677 | **672** |
| ctrl | 67 | 67 | 65 | **64** | **64** | 66 | 66 | 64 | 64 | **63** |
| router | 110 | 110 | 99 | **94** | 95 | 116 | 116 | 95 | **94** | 95 |

（DN 未列出，數字見上表；DN 在 32 組中只有 1 組最佳。）

32 組兩兩比較（贏:平:輸）：

| 比較 | 結果 |
|---|---|
| DLU 對 DI | 18:5:9 |
| DLA 對 DI | 18:4:10 |
| DLU 對 D / DR / DN | 26:2:4 / 26:2:4 / 28:2:2 |
| DLA 對 D / DR | 26:2:4 / 25:2:5 |
| DLA 對 DLU | 14:5:13 |

六種設定中最佳（平手都算）：DLA 14、DLU 10、DI 8、D 5、DR 5、DN 1。DLU 或 DLA 嚴格優於其他四種的有 18/32 組。

相對 D 的 AND 幾何平均比值：

| | DR | DI | DN | DLU | DLA |
|---|---|---|---|---|---|
| 全部 | 0.991 | **0.956** | 1.003 | 0.966 | 0.969 |
| 不含 arbiter | 0.994 | 0.982 | 1.000 | 0.971 | **0.970** |

- 局部 R 讓交替流程更好：勝場上 DLU / DLA 對 DI 約 18:9。DI 在幾何平均上領先，完全是 arbiter S1 的 785（lev 167）造成；去掉 arbiter 後 DLU / DLA 領先。
- 局部 R 版本明顯較好：bar（2720 / 2745，DI 2946 / 2952）、voter（7858 / 7750，DI 7928 / 8025）、router、ctrl、i2c、priority、square。
- 全電路 R 的 DI 較好：sin、sqrt，與 Phase 3d 一致。
- 仍由只跑 `&deepsyn` 的設定（D / DR）最好：cavlc、mem_ctrl；multiplier S2、int2float S1 也是 DR 最佳。
- level：局部 R 版本沒有出現 DI arbiter 那樣的 level 暴增；priority 的 level 反而大幅下降（85 → 37 / 39，99 → 36 / 35），max 也下降（約 330～420 → 276～287）。
- U 和 A 抽樣之間差距不大（DLA 對 DLU 14:13）。

### Phase 3d：局部 R 的正式比較（`phase3a.py --configs …`，`results/phase3d/`、`results/phase3d_deepsyn/`）

每個 job 1200 s，全部加 `-c`（只對被接受的結果做 cec），112 個 job 全部 cec ok。

**原始 EPFL 起點，seed 1**（3a 欄為同 seed 的全電路 R，`&deepsyn` 為 3a 的 20 分鐘結果）：

| 電路 | 3a U-e | 3a A-e | U-e-c | U-L2 | U-L4 | U-L8 | A-L4 | `&deepsyn` | iter 數（U-e-c / U-L4） |
|---|---|---|---|---|---|---|---|---|---|
| bar | 2949 | 2917 | 2951 | 2852 | 2831 | 2883 | **2768** | 2952 | 1966 / 9199 |
| log2 | 29006 | 28968 | 29006 | 29053 | 28892 | 28944 | 28881 | **28583** | 115 / 3976 |
| multiplier | 24301 | 24307 | 24301 | 24286 | 24290 | 24293 | **24277** | 24321 | 183 / 2812 |
| sqrt | 18249 | **18244** | 18249 | 18253 | 18251 | 18249 | 18252 | 18263 | 239 / 76438 |
| square | 15795 | 15793 | 15795 | 15793 | 15794 | 15791 | **15786** | 15798 | 280 / 10177 |
| voter | 7936 | 7924 | 7936 | 7936 | 7923 | **7756** | 7792 | 7935 | 628 / 42204 |
| sin | 4899 | 4896 | 4896 | 4891 | 4934 | 4907 | 4942 | **4888** | 605 / 25234 |
| max | 2779 | 2803 | 2779 | 2726 | 2724 | 2706 | 2724 | **2677** | 3917 / 13452 |
| priority | 379 | 367 | 377 | **365** | 366 | 366 | 368 | 376 | 17603 / 35302 |
| int2float | 151 | 144 | 151 | 165 | 145 | 143 | 144 | **122** | 17616 / 30126 |
| cavlc | 526 | 533 | 524 | 531 | 520 | 523 | 499 | **375** | 4748 / 11696 |
| arbiter | 8790 | 9869 | 8141 | 3646 | 3192 | 8655 | 3174 | **1998** | 604 / 28318 |
| mem_ctrl | 36487 | 36914 | 36492 | 39192 | 36285 | 35712 | 35700 | **7582** | 128 / 9472 |
| i2c | 847 | 828 | 846 | 842 | 789 | 809 | 778 | **684** | 7203 / 23562 |
| ctrl | 69 | 68 | 69 | 77 | 76 | **66** | 68 | 68 | 23899 / 17985 |
| router | 97 | 107 | 97 | 96 | 97 | 96 | **95** | 111 | 33970 / 43777 |

- A-L4 對全電路 R（U-e-c）贏 14/16，U-L4 贏 12/16；A-L4 對 3a A-e 贏 11/16。
- 大幅改善：bar 2951 → 2768（A-L4，也贏過 3a 所有設定的最佳 2911）、arbiter 8141 → 3174、voter 7936 → 7756（U-L8）、max 2779 → 2706、i2c 846 → 778、cavlc 524 → 499。
- 沒有幫助：sqrt、sin（需要全電路 R 的連鎖整理，次數多不能取代）。
- 對 `&deepsyn` 的勝負沒有翻盤：A-L4 贏 7/16（bar、multiplier、square、voter、sqrt、priority、router），log2、max、sin 和 3a 的「`&deepsyn` 大幅領先」電路仍然輸。
- h 沒有一致最佳值；A-L4 最穩定，L8 在 voter、max、ctrl 最好，L2 在 int2float、ctrl 明顯較差。
- level 大致不變（A-L4 對 U-e-c：arbiter 130 → 44、priority 79 → 59、max 183 → 227、sqrt 5982 → 6043）。

**`&deepsyn` 起點（3a seed 1 的 `&deepsyn` 結果），seed 2，與 3b（同起點、同 seed、全電路 R）比較**：

| 電路 | 起點 | 3b U-e | 3b A-e | 3b ILS 最佳（7 設定） | 3d U-L4 | 3d A-L4 | `&deepsyn` 繼續跑 |
|---|---|---|---|---|---|---|---|
| bar | 2952 | 2952 | 2911 | 2911 | 2800 | **2771** | 2952 |
| log2 | 28583 | 28443 | 28470 | 28433 | **28361** | 28392 | 28583 |
| multiplier | 24321 | 24309 | 24313 | 24298 | **24290** | 24296 | 24321 |
| sqrt | 18263 | **18253** | 18256 | 18253 | 18260 | 18260 | 18263 |
| square | 15798 | 15791 | 15791 | 15788 | **15783** | 15784 | 15798 |
| voter | 7935 | 7931 | 7927 | 7923 | 7919 | **7898** | 7927 |
| sin | 4888 | 4803 | 4800 | **4798** | 4826 | 4817 | 4888 |
| max | 2677 | 2676 | 2676 | 2676 | 2676 | 2676 | 2676 |
| priority | 376 | 370 | 371 | 370 | **368** | **368** | 376 |
| int2float | 122 | 119 | 117 | 117 | **114** | 117 | 122 |
| cavlc | 375 | **368** | 371 | **368** | 371 | 372 | 375 |
| arbiter | 1998 | 1975 | 1977 | 1975 | 1946 | 1953 | **1887** |
| mem_ctrl | 7582 | 7548 | 7553 | 7548 | 7509 | 7502 | **7370** |
| i2c | 684 | 681 | 678 | 678 | **662** | 666 | 675 |
| ctrl | 68 | 63 | 63 | 63 | 63 | 63 | 67 |
| router | 111 | 99 | 102 | 99 | 95 | **94** | 111 |

- 單一設定比較：U-L4 對 3b U-e 贏 11/16，A-L4 對 3b A-e 贏 10/16；U-L4、A-L4 對「`&deepsyn` 繼續跑」各贏 13/16。
- i2c 翻盤：3b 時 `&deepsyn` 繼續跑較好（675 對 678），局部 R 降到 662。
- 仍輸給 `&deepsyn` 繼續跑：arbiter、mem_ctrl。
- 全電路 R 仍較好：sin、sqrt、cavlc。

### Phase 2f：code 線在 R 之後的存活，以及強制保留（`phase2f.py`，`results/phase2f/`）

設定同 Phase 2e：原始 EPFL、R = rs3、K = 8、每組 200 次單次 kick、seed 1；每次 kick 跑兩種 R，kick 完全相同（`same_kicks` 全部 True）：
- **normal**：一般 R。
- **keep**（`-b`）：code 節點切開成 CO/CI，R 分別優化 E 與 D，再接回；r 條 code 線一定保留。

code 線是否存活以 2048 位模擬簽章判斷（不分正反相）；「新的」code 線指 x0 裡原本沒有的函數。

**存活比例（normal，16 個電路合計，不含 adder、dec）**：

| 抽樣方式 | code 線存活 | 新的 code 線存活 | 新 code 線全部消失的 kick |
|---|---|---|---|
| U | 36% | 32% | 67% |
| A | 50% | 36% | 67% |
| W4 | 47% | 42% | 60% |

**改善的 kick 幾乎都是 code 線被 R 消掉的那些**（normal，各抽樣方式合計，P(d_R < 0)）：

| 電路 | 新 code 線全部消失 | 至少一條新 code 線存活 |
|---|---|---|
| bar | 0.75（n=515） | 0.39（n=85） |
| ctrl | 0.86（133） | 0.18（465） |
| cavlc | 0.45（339） | 0.12（252） |
| mem_ctrl | 0.44（319） | 0.24（221） |
| priority | 0.33（417） | 0.01（176） |
| voter | 0.27（340） | 0.01（238） |
| arbiter | 0.25（558） | 0.00（41） |
| i2c | 0.25（455） | 0.03（144） |
| log2 | 0.26（379） | 0.06（213） |
| sqrt | 0.25（481） | 0.05（118） |

（其他電路方向相同。）這和 Phase 2e 的「只要一條新 code 線存活，92% 變差」一致。注意可能的混淆：容易被 R 完全吸收的 kick 本身可能比較小或比較便宜。

**強制保留（keep 減 normal，同一批 kick 的 d_R 中位數）**：
- 除了 dec W4（0），每個電路、每種抽樣方式（U、A、W4）都是正值，也就是保留 code 線讓結果變差：sqrt +17～21、square +21～43.5、priority +25～32、cavlc +29～34、ctrl +30～32、router +15.5～19。
- keep 的改善比例幾乎都是 0；例外是 bar（U：normal 76%、keep 75%，中位數 −58 對 −53）、mem_ctrl（keep 較好的比例約 30%，屬順序雜訊）、voter（keep 較好 23～30%，但 keep 的改善比例本身接近 0）。
- 結論：改善來自 R 把 E∘D 合併後重新整理，不是來自 r-bit 編碼本身是一個好的中間表示；「保留 r 條 code 線」不是有用的限制。bar 是唯一例外，它的編碼在保留時仍然有用。

## 6. Phase 3：Pareto search（go/no-go #2）

### 演算法

採 ε-constraint：對每個 delay 目標 D_k（取自 front B 的 delay 值），在 level ≤ D_k 下最小化面積。

```
x ← front B 中 level ≤ D_k 的最小面積 AIG
repeat until 預算用完:
    挑 window：以 slack(W) 加權（非 critical 區域權重高），critical 上的 window 只允許 min-level 候選
    挑抽樣方式與 code（依 Phase 2 的成功率做 bandit 式加權）
    y ← R(kick(x))
    若 level(y) ≤ D_k 且 area(y) < area(x)：x ← y
    否則以 SA 機率 exp(−Δ/τ) 接受（τ 逐步下降）
    每次更新 x 都更新 Pareto archive（含所有 D_k 共用）
```

- 每次接受的 AIG 都要先通過局部 `checkTruths`；寫入 archive 前跑全域 `cec`，不等價就中止並回報。
- Ablation：(a) 只用 local move 不 kick；(b) P0 取代重編碼；(c) 不用 SA（純 greedy）。
- 每組設定用 3 個 seed。

### 評估

- 每電路計算 archive front 與 baseline front B 的 **hypervolume 差**（M1、M2 各一）；參考點取兩者 area/delay 最大值 × 1.1。
- 在每個 D_k 回報 area 相對 B 的變化百分比。

**Go/no-go #2**：相同 wall-clock 下，若 hypervolume 沒有在 ≥ 3 個電路上超過 `&deepsyn` 加 3σ，結論就是
「本方法不優於 `&deepsyn`」，只做 Phase 4 的串接實驗後就結束。

**完成條件**：`phase3_pareto.csv`（每點：circuit, D_k, config, seed, time, ands, levels, map_area, map_delay, cec）、
`phase3_summary.csv`（每 circuit × config：hypervolume_M1, hypervolume_M2, Δarea@D_k）、所有 front 點的 AIG。
附 H3–H5 的結論。

## 7. Phase 4：與 `&deepsyn` 串接

以 `&deepsyn -T t` 的結果為起點，再以相同時間跑 Phase 3，與 `&deepsyn -T 2t` 比較。回答「本方法能不能當成
`&deepsyn` 卡住之後的補充擾動」。

## 8. 實作

新指令 `lsv_reencode`，作用在目前的 GIA 上（新檔 `src/ext-lsv/lsvReencode.cpp`，加進 `module.make`）：

```
lsv_reencode [-K num] [-M U|A|W|X|P0] [-s steps] [-N tries] [-I iters] [-D delay] [-T sec] [-S seed] [-R script] [-a] [-v]
  -a  只做 Phase 1 統計，不修改電路
```

- 從 `lsvPipeAcd.cpp` 抽出可共用的部分（window、image、合成、cec）放到共用的 header/檔案，不要複製貼上。
- R 透過 `Cmd_CommandExecute` 呼叫既有 ABC 指令，不重寫最佳化。
- window 輸出邊界、隨機 code 產生器（U、A、W、X、P0，含正規化去重）、SA 迴圈、archive 是新寫的部分。
- 隨機數一律由 `-S` seed 決定，相同 seed 必須產生相同結果（runtime 預算改用固定迭代數時）。
- 實驗 driver 用 Python script 放在 `src/develop/window_reencode/`，可一鍵重跑全部 CSV。

### 單元測試（先寫，再實作）

| 測試 | 預期 |
|---|---|
| 3-to-8 decoder，window = 整個電路 | N = 8，r = 3，`|G| > r` |
| 對每種抽樣方式隨機 1000 個 code，在小電路上套用 | 全部 `checkTruths` 與 `cec` 通過 |
| G 刻意不取輸出邊界 | 被偵測並拒絕（舊邏輯無法刪除） |
| 同一 seed 跑兩次 | AIG 完全相同 |
| 刻意塞入錯誤 decoder | `checkTruths` 擋下；跳過 `checkTruths` 時 `cec` 擋下，流程中止 |
| N = 5, 6, 8 時 X 窮舉出的編碼數 | 分別為 140、420、840 |
| 只做 code 線排列/反相的 φ | 正規化後與原 code 相同（被去重） |

## 9. 範圍外（這次不做）

- 冗餘 code（m > r 條線）。
- 使用下游 observability don't care 的 class 合併（column multiplicity）。若 Phase 2 顯示 image-based class 的擾動
  空間太小，再另外評估。
- 跨 window 串接 code、不 decode 回 G 的全域重編碼。
- Sequential 電路、retiming。
- hyp：cec 在 3600 秒內跑不完，無法驗證（Phase 0 已確認）。
- 修改 ABC 內建最佳化指令本身。

## 10. 執行注意事項

- 先在 ctrl、int2float、cavlc、router、priority 上跑通，再擴大到 arithmetic 電路。hyp、div、sqrt 很大，`dc2` 一次
  就要數秒，迭代數會很少，結果要另外標註。
- 比較一律在**相同 wall-clock** 下，不比迭代數。
- 單一電路超過預算 2 倍就中止並回報，不要默默縮小 K 或迭代數。
- 任何改善小於 max(σ, 0.5%) 的結果不得宣稱為改善（B 群多數電路 σ = 0）。
- 與本計劃假設不符的發現，優先回報。
