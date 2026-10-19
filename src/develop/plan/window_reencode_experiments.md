# Window 重編碼實驗：設計總表

這份文件整理目前所有實驗的設計：baseline、執行指令、預算、方法與結果，以及每項比較的限制。
- 結果細節與討論在 `window_reencode_pareto.md`。
- 程式：`src/ext-lsv/lsvReencode.cpp`（`lsv_reencode`）、`src/develop/window_reencode/phase*.py`。
- 原始數據：`src/develop/window_reencode/results/`。

## 1. 共通設定

### 環境與 benchmark

| 項目 | 設定 |
|---|---|
| 電路 | EPFL combinational suite，`../testcases/aig/epfl/<c>.aig`（binary AIGER）。Phase 2c 之後使用 18 個電路：adder、bar、log2、multiplier、sqrt、square、voter、sin、max、priority、int2float、cavlc、arbiter、mem_ctrl、i2c、ctrl、router、dec；Phase 3a 的多 seed 實驗與 3b、3c 排除 adder、dec（在所有方法下都不變） |
| 工具 | 本 repo 的 `./abc`（branch `feature-nf-flag-Y`，含未 commit 的 `lsv_reencode`） |
| 機器 | 42 核心的共用伺服器（有其他使用者）。每個 job 是單執行緒的 `abc` 程序，40～43 個 job 並行 |
| 時間 | `lsv_reencode -T` 與 `&deepsyn -T` 都是 wall clock（`Abc_Clock`，CLOCK_MONOTONIC），每輪結束後才檢查，所以最多超時一輪 |

### 指標

| 指標 | 定義 |
|---|---|
| 面積 | AIG 的 AND 節點數（`Gia_ManAndNum`），**不是** mapping 後的面積。asap7 mapping（M2）只在 Phase 0 做過 |
| delay | AIG level（unit delay）；目前只記錄，不設限制 |
| d_R（Phase 2 系列） | AND(R(kick(x))) − AND(R(x))；< 0 代表 kick 加上 R 之後比 R 直接作用在 x 上更小 |
| 正確性 | 每次 kick 都對 x 做 cec（`cecGiaSafe`：在 fork 的子程序跑 Fraig CEC，失敗則退回 SAT sweeping），回傳 1（等價）、0（不等價）或 −1（undecided）。每個 window 也會窮舉檢查 E∘D = F_W（`composesTo`）。ILS 每次 iteration 對目前電路做 cec，結束時再對 x0 做一次；最終結果另外用 `&cec` 獨立確認 |

### Re-optimization script R

| 名稱 | 內容 | 使用階段 |
|---|---|---|
| dc2 | `&put; dc2; &get -n` | Phase 2 |
| dc2x5 | `&put; dc2` ×5 `; &get -n` | Phase 2b |
| **rs3**（主要） | `&put; C; C; C; &get -n`。C 為 `compress2rs`：`balance -l; resub -K 6 -l; rewrite -l; resub -K 6 -N 2 -l; refactor -l; resub -K 8 -l; balance -l; resub -K 8 -N 2 -l; rewrite -l; resub -K 10 -l; rewrite -z -l; resub -K 10 -N 2 -l; balance -l; resub -K 12 -l; refactor -z -l; resub -K 12 -N 2 -l; rewrite -z -l; balance -l`，與 `&deepsyn` 內部使用的字串相同 | Phase 2b 起 |

### Window

- **root**：從可到達 PO 的 AND 節點中均勻抽一個。
- **L**：從 root 往上做 reconvergence-driven 展開（`findWindow`），|L| ≤ K。
- **C**：L 的 forward closure，上限 500 個 AND（`-C`）。
- **G**：C 中被 C 外部或 PO 使用的節點（輸出邊界）。
- **cyclic 規則**：被 leaf 依賴的 C 節點，以及 C 內依賴它們的節點，都從 C 移除；若 root 被移除，這個 window 作廢。
- **class**：依 G 的值把 2^|L| 個 leaf pattern 分組，共 N 組，r = ⌈log₂N⌉。只接受 N ≥ 5 的 window（`-N 5`）。
- **重試**：每次 kick 最多嘗試 200 個 root，直到找到合格且合成成功的 window。

### 抽樣方式（決定 class → codeword 的指派 φ）

| 代號 | 方式 |
|---|---|
| U | 均勻隨機 φ |
| A | φ = M·base ⊕ c；M 是隨機可逆、但不是排列的矩陣；base 是 `assignCode` 的貪婪編碼 |
| W1 / W4 | 從 base 出發，做 1 次或 4 次「交換兩個 class 的 codeword」或「移到未使用的 codeword」 |
| base | `assignCode` 的確定性編碼（對照組） |
| P0 | 不經 code 瓶頸，直接從真值表合成 L → G（對照組，結果確定） |
| N | null kick：把 C 原封不動放回，只改變節點編號順序（對照組，用來量 R 的順序雜訊） |
| X | N ≤ 8 時窮舉所有實質不同的編碼（只用於 Phase 1） |

φ 會先正規化（class 0 → 0、code 欄排序），再進行合成。

### 合成與預算

- **E 的合成**：真值表是完全指定的，所以只有 1 種填法。
- **D 的合成**：沒用到的 codeword 是 don't care，有 3 種填法：填 0、填 1、ISOP。
- **候選**：每種填法各做 Kit 因式分解（`Kit_TruthToGia`）、`Gia_ManCompress2`（dc2）、`Gia_ManAigSyn2`（syn2）三種，取 (AND, level) 最小者。
- **預算**：每個電路（E、D 或 P0 的 L → G）各自的上限是 B·(|C| + 10)，預設 B = 3（`-B`；0 表示不設上限）。超過上限就放棄這個 code，換下一個 root。檢查點是 Kit 因式分解的原始結果。
- **替換**：重建 AIG 的順序是：先放 window 下游以外的節點，再放 E/D，最後重建下游；過程中開啟 structural hashing。所以即使 window 內容不變，節點編號也會改變，這就是 N 量到的順序效應。

### 執行指令範本

```
# 單次 kick（Phase 2 系列）：從同一個 x 出發，每次獨立
&r <x>; lsv_reencode -k -M <sampler> -K 8 -n 200 -S 1 -R "<R>" -o <out.csv> [-m 1 | -P 10] [-B 3] [-W <dir>] [-w <best.aig>]
# ILS（Phase 3）：接受改善後從新電路繼續
&r <x>; lsv_reencode -i -M <sampler> [-e] -K 8 -n 1000000 -T <sec> -S <seed> -R "<rs3>" -o <out.csv> -w <final.aig>
# 對照組
&r <x>; &deepsyn -T <sec> -S <seed> -v; &w <out.aig>
```

ILS 的流程：
1. cur = min(x, R(x))，以 (AND, level) 字典序比較。
2. 每次 iteration：kick(cur) → cec → R → 得到 z。
3. 接受條件：AND(z) < AND(cur)，或 AND 相同但 level 較低；加上 `-e` 時，AND 相同且 level 不變差也接受（平移）。

`&deepsyn` 使用預設值：`-I 1`、`-J 10⁹`、不加 `-c`、`-t`、`-o`。每一輪的內容是：
1. `&dch [-f]`
2. `&if -a -K k`（k 在 3～6 之間輪替）
3. `&mfs -e -W 20 -L 20`
4. 隨機選擇是否做 `&fx; &st`
5. `&dc2` 或 `compress2rs` 擇一

它以 random walk 前進，只記錄看過的最佳解。

## 2. 各階段

| 階段 | 問題 | 起點 | 方法 / 對照 | 規模 | 預算 | Driver |
|---|---|---|---|---|---|---|
| 0 | baseline 與雜訊 | 原始 EPFL（20 個） | dc2、resyn2、compress2rs（各自反覆到收斂）、balance、`&syn2`、`&b -d`、`&deepsyn -T {60, 600} -S {1, 2, 3}`；asap7 `&nf -p [-R 0..100]` | 每電路約 13 個變體 | `&deepsyn` 60 s / 600 s | `phase0.py` |
| 1 | 直接 move 的 gain（H1） | Phase 0 中最小的 AIG（12 個） | `lsv_reencode -a`：P0、base、U、A、W1、W4 各 16 個 code，X | 最多 1000 個 root，K ∈ {6, 8, 10, 12} | B = 3 | `phase1.py` |
| 2 | 單次 kick + R（go/no-go #1） | 同 Phase 1 | U、A、W1、W4，對照 base、P0；R = dc2 | 200 次 kick × 12 電路 × 6 種 × K ∈ {6, 8} | — | `phase2.py` |
| 2b | R 太弱嗎？ | 同上（log2、sin、sqrt、priority） | 同上；R = dc2x5、rs3 | 200 次 × K = 8 | — | `phase2b.py` |
| 2c | 從未最佳化的電路出發，結構能變多大 | 原始 EPFL（18 個） | U、A、W4、base、P0；R = rs3；m1（1 個 window × 200 次）、P10 / P50（替換到 10% / 50% 的 AND × 20 次） | 21,446 次 kick | B = 3 | `phase2c.py` |
| 2d | 放大 window | 原始 EPFL（8 個） | U、P0，K ∈ {6, 8, 10, 12} | 200 次 | B ∈ {3, 10, 0} | `phase2d.py` |
| 2e | 順序雜訊與殘留分析 | 原始 EPFL（18 個） | N、P0、base、U、A、W4；R = rs3；`-W` 追蹤節點來源 + `residue.py` | 21,476 次 kick | B = 3 | `phase2e.py` |
| 3a | ILS 對 `&deepsyn` | 原始 EPFL | U、A、P0、N 及其 -e 變體；對照組 `&deepsyn` | 18 個電路 × seed 1；16 個電路 × seed 2、3（不含 N） | 每個 job 1200 s | `phase3a.py` |
| 3b | ILS 疊在 `&deepsyn` 之後 | 3a 中 `&deepsyn` seed 1 的結果 | 同 3a，seed 2；對照組 `&deepsyn -S 2` 從同一起點繼續跑 | 16 個電路 | 每個 job 1200 s | `phase3a.py --start` |
| 3c | 交替流程（完成） | 原始 EPFL | D、DR、DI、DN（見第 3 節） | 16 個電路 × seed 1、2 | 每個 job 3600 s | `phase3c.py` |
| 3d | 局部 R（完成） | 原始 EPFL；另加 3a 中 `&deepsyn` seed 1 的結果 | U-e 全電路 R 對照、U-e 局部 R（`-L` 2/4/8）、A-e 局部 R（`-L` 4）；全部加 `-c`（只對被接受的結果做 cec） | 16 個電路，seed 1（`&deepsyn` 起點用 seed 2） | 每個 job 1200 s | `phase3a.py --configs …` |
| 3c 追加 | 交替流程改用局部 R（完成） | 原始 EPFL | DLU：(`&deepsyn -T 600`; ILS `-M U -e -c -L 4`) × 3；DLA：同上但 `-M A` | 16 個電路 × seed 1、2 | 每個 job 3600 s | `phase3c.py --configs DLU DLA` |
| 2f | code 線的存活與保留（完成） | 原始 EPFL | U、A、W4、base，每次 kick 跑兩種 R：一般 R，以及 `-b`（code 線切開成 PO/PI，R 不能跨過它，之後再接回）；兩種模式的 kick 完全相同，可以逐次配對 | 18 個電路 × 200 次 | — | `phase2f.py`（使用 `./abc_dev`） |

### 各階段結果摘要

- **Phase 0**：`&deepsyn -T 600`（3 個 seed 取最小）相對於 compress2rs 反覆到收斂的結果：
  - 大幅領先：arbiter −84%、mem_ctrl −80%、int2float −42%、cavlc −39%、i2c −33%、ctrl −27%、router −23%、priority −13%。
  - 幾乎卡住：bar −6%、max −5%、log2 −1.5%、sin −1.6%，adder、multiplier、sqrt、square、voter 介於 0 和 −0.1%。
  - 不同 seed 間的差距記為 σ，作為 Phase 2 門檻 max(σ, 0.5%) 的一部分。
- **Phase 1**：H1 成立。
  - 任何抽樣方式的直接 gain > 0 比例都 < 1%（例如 K = 8 時 P0 0.16%、U 0.03%）。
  - U 的合成失敗率：K = 6 時 25%，K = 10、12 時 67～69%。
- **Phase 2**：嚴格判定 NO-GO。
  - undo 率 < 80%，但最佳改善都低於門檻（log2 −6、sin −3…）。
  - 對照組 P0 不輸重編碼。
  - Phase 2e 證明 log2 和 sin 的改善落在節點順序雜訊的範圍內。
- **Phase 2b**：rs3 比較能吸收 kick（undo 38～48%），但增益仍然 ≤ 0.15%，對照組也不輸。
- **Phase 2c**：
  - 結構可以大幅改變（P50 的 new_frac 中位數約 0.44），但面積 +37～59%。
  - m1 中確實超越 R* 的改善：bar −177、priority −14、router −14、ctrl −11。
- **Phase 2d**：大 window 的隨機編碼多半超出預算；放寬預算後電路膨脹，R 救不回來。最佳結果都來自 |L| = 3～4 的 window。
- **Phase 2e**：
  - null kick 的雜訊因電路而異：9 個電路為 0，mem_ctrl 在 −94～+91 之間。
  - 只要有一條「新的」code 線在 R 之後存活，92% 變差。
  - 「XOR 型 code 難消」是不同電路混在一起造成的混淆，已更正。
- **Phase 3a / 3b**：見第 4 節。
- **Phase 3c**（128 個 job 全部 cec ok）：DI 對 D 24:5、對 DR 23:6、對 DN 25:5（其餘平手），32 組中 22 組最佳；DN 對 DR 11:12。增益來自重新編碼的 kick，不是重開或順序雜訊。DI 輸在 mem_ctrl、cavlc、int2float。arbiter S1 的 DI 785 level 升到 167。
- **Phase 3d**：局部 R 讓 iteration 數增加 5～300 倍。原始起點下 A-L4 對全電路 R 贏 14/16（bar 2951 → 2768、arbiter 8141 → 3174、voter 7756（L8）、i2c 778）；sqrt、sin 沒有幫助；對 `&deepsyn` 仍是 7/16。`&deepsyn` 起點下 U-L4 / A-L4 對「`&deepsyn` 繼續跑」各贏 13/16（i2c 翻盤為 662），仍輸 arbiter、mem_ctrl。
- **Phase 3c 追加**（64 個 job 全部 cec ok）：DLU / DLA 對 DI 18:9 / 18:10，對 D 都是 26:4；六種設定中 DLA 14 組最佳、DLU 10、DI 8。幾何平均（相對 D）DI 0.956 最好，但只因為 arbiter S1 的 785（lev 167）；不含 arbiter 時 DLA 0.970、DLU 0.971、DI 0.982。sin、sqrt 仍是 DI 較好，cavlc、mem_ctrl 仍是 D / DR 較好。
- **Phase 2f**：一般 R 之後新 code 線只有 32～42% 存活，約 2/3 的 kick 中全部消失；改善幾乎都發生在新 code 線全部消失的 kick（例如 priority 0.33 對 0.01）。強制保留（`-b`）在所有電路上中位數都變差（+4.5～+43.5，dec 約 0），只有 bar 仍保有改善。保留 r 條 code 線不是有用的限制。

## 3. Phase 3c 的設計

| 代號 | 流程（總長 3600 s，每一段從上一段的結果繼續） | 目的 |
|---|---|---|
| D | `&deepsyn -T 3600 -S 10s` | 基準 |
| DR | `&deepsyn -T 600` × 6，seed 依序為 10s + k | 只重啟 `&deepsyn` 的效果 |
| DI | (`&deepsyn -T 600`; ILS `-M U -e -T 600`) × 3 | 交替流程 |
| DN | (`&deepsyn -T 600`; ILS `-M N -e -T 600`) × 3 | 對照組：只有 R 與順序雜訊，沒有 kick |
| DLU | (`&deepsyn -T 600`; ILS `-M U -e -c -L 4 -T 600`) × 3 | 追加：DI 改用局部 R |
| DLA | 同 DLU，但 `-M A` | 追加：Phase 3d 最穩定的局部 R 設定 |

- s 為 seed（1、2），k 為段的編號。四種設定的第 1 段完全相同（同為 `&deepsyn -T 600 -S 10s`）。
- ILS 的 seed 為 100s + k；K = 8，R = rs3。
- 輸出在 `results/phase3c/`：
  - `steps.csv`：每段的輸入/輸出 AND 與 level、秒數。
  - `phase3c_summary.csv`：最終值，以及每 600 秒邊界時的 AND 數。
- 每個 job 的最終結果都用 `&cec` 和 x0 比對。
- 判讀方式：
  - DI < D 且 DI < DR：交替有效，而且不只是重啟的效果。
  - DI < DN：增益來自 kick，不是 R 或順序雜訊。
- 已知限制：
  - DI 固定使用 U-e。這是根據 3a 的幾何平均比值選的（1.264，最好），但若以勝場數計，A-e 較好（3a 48 場中贏 21 場，U-e 只贏 15 場）。
  - 只有 2 個 seed。
  - 啟動時機器上已有其他使用者的負載（load 約 14），加上 43 個並行 job，會有 CPU 競爭。

## 4. Phase 3a / 3b 的公平比較（更正）

`window_reencode_pareto.md` 中 Phase 3a/3b 表格的「ILS 最佳」是 6～8 種 ILS 設定（各 1200 s）取最小值，`&deepsyn` 則只跑一次 1200 s。所以 ILS 等於用了 6～8 倍的運算量，**表格中的勝幅因此偏大**。

改成每種設定單獨和同 seed 的 `&deepsyn` 比較：

**Phase 3a**（原始起點，seed 1～3，每格是 3 個 seed 中贏 `&deepsyn` 的次數）：

| 電路 | U | A | P0 | U-e | A-e | P0-e |
|---|---|---|---|---|---|---|
| router | 2 | 2 | 0 | 3 | 3 | 0 |
| priority | 3 | 3 | 0 | 2 | 3 | 0 |
| bar | 0 | 2 | 0 | 1 | 3 | 0 |
| voter | 2 | 3 | 3 | 1 | 3 | 3 |
| sqrt | 2 | 2 | 3 | 2 | 2 | 3 |
| sin | 1 | 2 | 2 | 1 | 2 | 2 |
| multiplier | 3 | 3 | 3 | 3 | 3 | 3 |
| square | 2 | 2 | 2 | 2 | 2 | 2 |
| ctrl | 0 | 1 | 0 | 0 | 0 | 0 |
| log2、max、i2c、int2float、cavlc、arbiter、mem_ctrl | 0 | 0 | 0 | 0 | 0 | 0 |
| 合計（48 場） | 15 | 20 | 13 | 15 | **21** | 13 |
| 至少贏 2 個 seed 的電路數 | 6 | **8** | 5 | 5 | **8** | 5 |

**Phase 3b**（從 `&deepsyn` 的結果出發，單一 seed，和繼續跑 20 分鐘的 `&deepsyn` 比較）：

| 設定 | 贏 | 平手 | 贏的電路 |
|---|---|---|---|
| A、A-e | 11 | 1 / 2 | router、priority、ctrl、bar、sqrt、sin、multiplier、square、log2、int2float、cavlc |
| U、U-e | 10 | 2 | 同上，但不含 bar |
| P0 | 8 | 3 | priority、bar、voter、sqrt、sin、multiplier、square、log2 |
| N（null kick） | 4 | 6 | sin、multiplier、square、log2 |

結論：
- **單一固定設定**（A 或 A-e）仍然贏過 `&deepsyn`：原始起點下有 8 個電路，`&deepsyn` 之後有 11 個電路。
- N 也贏了 sin、multiplier、square、log2，代表這 4 個電路的增益有一部分是順序雜訊。
- router、priority、ctrl、bar、int2float、cavlc 的增益則來自 kick（N 在這些電路上沒有任何改善）。

## 5. 效度上的限制（判讀結果時要記住）

1. **面積指標是 AND 數**，不是 mapping 後的面積；delay 沒有限制，部分結果 level 會升高（例如 3b 的 int2float 13 → 16）。
2. **seed 數很少**：3a 有 3 個 seed，3b、3d 只有 1 個，3c 有 2 個。勝負在 1～2 個 AND 以內的電路（square、max、ctrl）都不能下定論。
3. **R 的節點順序雜訊**：在 mem_ctrl、log2、i2c、sin、cavlc 等電路上，改善要和 N 比較才有意義。
4. **wall clock 與 CPU 競爭**：機器是共用的，各設定交錯執行、受到相同的影響。但絕對時間不能拿來和其他論文的數字比較。
5. **`&deepsyn` 只用預設參數**，沒有調整 `-I`、`-J` 或 choices；Phase 0 的其他變體（例如 bar 的 `&syn2` 2688）在某些電路上比 `&deepsyn` 更好。
6. **ILS 每次 iteration 都對整個電路跑 R**（3d 之前），大電路在 20 分鐘內只有 100～200 次 iteration，結果偏向保守；3d 的局部 R（`-L`）解決了次數問題，但在 sqrt、sin 上反而較差。
