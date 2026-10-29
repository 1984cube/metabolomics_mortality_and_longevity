[← 上一篇：07 · 基因程序评分与跨组织比对](07-gene-program-scoring.md) · [目录](README.md) · [下一篇：09 · 个体级队列分析：CHIP × 骨密度 × 骨折，以及共病 →](09-ukb-cohort-and-comorbidity.md)

---

# 08 · 虚拟扰动：CellOracle、scTenifoldKnk、NicheNet 与 Geneformer

> 这一章带走什么：四种工具各自在做什么；为什么复杂模型必须和线性基线并排比较；"有没有方向正确的先验网络"为什么决定成败；怎么用已知答案回测。

### 9.1 要回答的问题

报告想检验的假说里有一个"可干预的节点"：抑制它，应该**同时**降低骨髓里的破骨样评分和斑块里的钙化评分。但谁来做这个实验？报告的第一步是**虚拟扰动**：用计算模型预测一个基因被敲除或抑制后，细胞状态会怎样变。

报告给出的工具分工：

- 转录因子：CellOracle；
- 非转录因子：scTenifoldKnk；
- 分泌因子：NicheNet 的配体活性；
- Geneformer：只做人类数据的推理，并与线性基线并排报告；
- 用已有答案做回测：TET2、MMP9、OSM、IL-20 都有已发表的敲除或中和实验。

### 9.2 四种工具各做什么

| 工具 | 输入 | 核心思路 | 擅长 | 局限 |
|---|---|---|---|---|
| CellOracle（Kamimoto 等，Nature 2023） | scRNA-seq，加上 scATAC-seq 或基序给出的"底图网络" | 先推断每种细胞类型自己的调控网络（转录因子到靶基因的系数），再把扰动沿网络一轮一轮传播，预测细胞状态变化 | 转录因子的敲除 | 依赖底图网络；线性假设；对间接调控和反馈有限 |
| scTenifoldKnk（Osorio 等，Patterns 2022） | 只需 scRNA-seq | 用单细胞数据构建基因网络，"敲除"某个基因就是把它发出的边设为 0，比较前后网络，找受影响的基因 | 没有转录因子注释的基因 | 网络本质是相关性网络，没有方向先验 |
| NicheNet（Browaeys 等，Nat Methods 2020） | 发送细胞、接收细胞的表达数据，加上"配体到靶基因"的先验知识 | 不预测敲除，而是回答"哪个配体最能解释接收细胞里观察到的基因变化"，输出配体活性排序 | 分泌因子、细胞间通讯 | 配体必须被发送细胞表达、受体必须被接收细胞表达；先验知识有偏 |
| Geneformer（Theodoris 等，Nature 2023） | 单细胞转录组 | 在约 3,000 万个人类单细胞转录组上预训练的 Transformer；把细胞变成按相对表达量排序的基因序列；"敲除"就是从序列里删掉这个基因，看细胞表征往哪里移动 | 无需自己建网络，可以一次测很多基因 | 黑箱；预测未必优于简单的线性基线 |

以上是概述，具体设计（回归方式、传播轮数、网络构建细节）以各自原文和文档为准。

### 9.3 实战 A：玩具基因调控网络和扰动传播

用一个"真相已知"的小网络演示原理。设计：

- 40 个基因，其中前 5 个是转录因子（TF）；
- TF1 调控 TF2（所以 TF2 不是"根"）；TF3 调控 TF5；
- 每个靶基因受 1 到 2 个 TF 直接调控；一部分靶基因还受其他靶基因调控（间接调控）；
- 真相：把某个基因敲到 0，其余基因的平均表达变化多少（由模型精确算出）。

```r
set.seed(61)
p <- 40; tf <- 1:5                                     # 40 个基因，其中前 5 个是转录因子（TF）
B <- matrix(0, p, p)                                   # B[j, k]：基因 k 对基因 j 的直接调控强度
B[2, 1] <- 0.9                                         # TF1 调控 TF2（TF2 不是“根”）
B[5, 3] <- 0.8                                         # TF3 调控 TF5
for (j in 6:p) {                                       # 每个靶基因受 1–2 个 TF 直接调控
  reg <- sample(tf, sample(1:2, 1))
  B[j, reg] <- sample(c(-1, 1), length(reg), TRUE) * runif(length(reg), 0.4, 0.9)
}
for (j in 21:p) B[j, sample(6:(j - 1), 1)] <- sample(c(-1, 1), 1) * runif(1, 0.3, 0.6)   # 靶基因之间的间接调控
mu <- c(rep(3, 5), rep(2, p - 5))                      # 每个基因的平均表达量
I_B <- diag(p) - B
cvec <- as.vector(I_B %*% mu)                          # 让平均表达量恰好等于 mu 的截距项
n <- 3000
X <- (matrix(cvec, n, p, byrow = TRUE) + matrix(rnorm(n * p), n, p)) %*% t(solve(I_B))   # x = (I-B)^-1 (c + e)
colnames(X) <- paste0("G", 1:p)

# 真相：把基因 k 敲除到 0（do 操作），其余基因的平均表达量变化多少
true_ko <- function(k) {
  Bk <- B; Bk[k, ] <- 0; ck <- cvec; ck[k] <- 0
  (solve(diag(p) - Bk) %*% ck)[, 1] - mu                # 敲除后的平均值 - 原来的平均值
}
```

三种预测方法：

- **线性基线：** 每个基因单独对被敲除的 TF 做回归，斜率乘以"敲到 0"的变化量；
- **传播（无先验）：** 任何 TF 都可以是任何基因的候选调控因子，用岭回归推断系数，再把扰动沿网络传播 3 轮（CellOracle 式思路的教学简化版，不是 CellOracle 本身）；
- **传播（有先验）：** 只允许"底图网络"里的边（相当于 scATAC-seq 给出的方向正确的候选边：真实的边都在，再混入 20% 的假候选边，但方向不会反）。

```r
# 做法 A（CellOracle 式思路的教学简化版）：先推断“调控因子 -> 基因”的系数，再把扰动沿网络传播 3 轮
ridge <- function(y, Z, lambda = 1) solve(crossprod(Z) + lambda * diag(ncol(Z)), crossprod(Z, y))
fit_grn <- function(cand) {                             # cand[j, k] = TRUE：允许基因 k 作为基因 j 的调控因子
  Bhat <- matrix(0, p, p)
  for (j in 1:p) if (any(cand[j, ]))
    Bhat[j, cand[j, ]] <- ridge(X[, j] - mean(X[, j]), scale(X[, cand[j, ], drop = FALSE], scale = FALSE))
  Bhat
}
cand_free <- matrix(FALSE, p, p); cand_free[, tf] <- TRUE; diag(cand_free) <- FALSE   # 无先验：任何 TF 都可能调控任何基因

cand_prior <- matrix(FALSE, p, p)                       # 有先验（相当于 scATAC-seq 给出的“底图网络”）：
cand_prior[, tf] <- (B != 0)[, tf]                      #   真实的 TF -> 基因 边都在
for (j in 1:p) for (k in tf) if (j != k && B[j, k] == 0 && B[k, j] == 0 && runif(1) < 0.2) cand_prior[j, k] <- TRUE   # 再混入 20% 的假候选边，但方向不会反

propagate <- function(Bhat, k, steps = 3) {
  delta <- numeric(p); delta[k] <- -mean(X[, k]); total <- delta
  for (s in seq_len(steps)) { delta <- as.vector(Bhat %*% delta); delta[k] <- 0; total <- total + delta }
  total
}
# 做法 B（线性基线）：每个基因单独对 TF k 做回归，斜率 × (0 - 平均表达量)
baseline <- function(k) sapply(1:p, function(j) cov(X[, j], X[, k]) / var(X[, k])) * (0 - mean(X[, k]))

auc <- function(score, truth) {                         # 用打分把“真正受影响的基因”排到前面的能力
  r <- rank(score); n1 <- sum(truth); n0 <- sum(!truth)
  (sum(r[truth]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
B_free <- fit_grn(cand_free); B_prior <- fit_grn(cand_prior)
evaluate <- function(k) {
  tr <- true_ko(k)[-k]; affected <- abs(tr) > 0.2       # 真相里变化超过 0.2 的基因
  pred <- list(线性基线 = baseline(k)[-k], 传播_无先验 = propagate(B_free, k)[-k], 传播_有先验 = propagate(B_prior, k)[-k])
  list(n = sum(affected), cor = sapply(pred, function(v) cor(tr, v)), auc = sapply(pred, function(v) auc(abs(v), affected)))
}
e1 <- evaluate(1); e2 <- evaluate(2)
cat("真相里受影响的基因数：TF1 =", e1$n, "；TF2 =", e2$n, "\n")
cat("\n预测与真相的相关系数\n");  round(rbind(`TF1（上游没有调控因子）` = e1$cor, `TF2（被 TF1 调控）` = e2$cor), 2)
cat("\n把受影响基因排到前面的能力（AUC）\n"); round(rbind(`TF1（上游没有调控因子）` = e1$auc, `TF2（被 TF1 调控）` = e2$auc), 2)
```

```text
真相里受影响的基因数：TF1 = 20 ；TF2 = 12 

预测与真相的相关系数
                        线性基线 传播_无先验 传播_有先验
TF1（上游没有调控因子）     1.00        1.00        0.99
TF2（被 TF1 调控）          0.72        0.73        0.97

把受影响基因排到前面的能力（AUC）
                        线性基线 传播_无先验 传播_有先验
TF1（上游没有调控因子）     1.00        1.00        0.96
TF2（被 TF1 调控）          0.82        0.82        0.88
```

读法：

1. **对没有上游调控因子的 TF1，线性基线已经几乎完美**（相关 1.00）。复杂模型没有用武之地。这就是为什么必须和线性基线并排比：它经常已经够好。
2. **对被 TF1 调控的 TF2，线性基线和无先验传播都只有 0.72 和 0.73。** 原因：TF1 的靶基因和 TF2 都相关（因为 TF2 本身受 TF1 驱动），两种方法都把这些基因误当成 TF2 的靶基因。调试时我还看到，无先验的回归会推断出"TF2 调控 TF1"这条反向边，于是扰动传到了 TF1 的靶基因上。
3. **有了方向正确的先验网络，相关升到 0.97，AUC 升到 0.88。** 这就是 CellOracle 要借助 scATAC-seq 给出"底图网络"的原因：它告诉模型哪些转录因子可能结合哪些基因，方向由结合位点固定。
4. 即使有先验，AUC 也只有 0.88：靶基因之间的间接调控没有被建模，仍有遗漏。

### 9.4 实战 B：NicheNet 的思路

NicheNet 不预测"敲除后会怎样"，而是反过来：**接收细胞里观察到了一批差异表达基因，哪个配体最能解释它们？**

玩具设计：30 个配体、300 个靶基因；先验知识是一个"配体到靶基因的调控潜力"矩阵；真正起作用的是配体 7，接收细胞里差异表达的基因是配体 7 的靶基因的 70%，再混入 15 个无关的噪声基因。每个配体的"活性"取它的先验潜力向量与"是否差异表达"之间的 Pearson 相关（NicheNet 默认用来给配体排序的度量）。

```r
set.seed(62)
n_lig <- 30; n_tar <- 300
P <- matrix(0, n_tar, n_lig)                            # 先验知识：配体 -> 靶基因的“调控潜力”
for (l in 1:n_lig) { tg <- sample(n_tar, 25); P[tg, l] <- runif(25, 0.2, 1) }
colnames(P) <- paste0("配体", 1:n_lig)
true_l <- 7                                             # 真正在起作用的配体
tg7 <- which(P[, true_l] > 0)
de <- integer(n_tar)                                    # 受体细胞里观察到的差异表达基因
de[sample(tg7, round(0.7 * length(tg7)))] <- 1L         # 配体 7 的靶基因里有 70% 被检出
de[sample(which(de == 0), 15)] <- 1L                    # 再混入 15 个无关的噪声基因
activity <- apply(P, 2, function(v) cor(v, de))         # NicheNet 的配体活性：先验潜力与差异基因的 Pearson 相关
round(sort(activity, decreasing = TRUE)[1:5], 3)
```

```text
 配体7 配体26  配体3 配体16 配体23 
 0.538  0.091  0.049  0.043  0.034 
```

配体 7 的活性是 0.538，排第一，第二名只有 0.091。

真实 NicheNet 的使用要点：

1. 先验模型整合了配体到受体、信号通路、基因调控三类数据库；
2. 配体必须在发送细胞里表达，受体必须在接收细胞里表达，否则"活性"没有意义；
3. 用于报告：发送细胞是 CHIP 突变单核细胞，接收细胞是瓣膜间质细胞，看 OSM 是否排在前面。报告说 OSM 的作用已有发表，这正好是一个"已知答案回测"。

### 9.5 实战 C：Geneformer 吃进去的是什么

Geneformer 的输入不是数字矩阵，而是**每个细胞一条按相对表达量排序的基因序列**。"相对"是指每个基因除以它在整个训练语料库里的典型表达水平，这样管家基因（GAPDH、ACTB）被压到末尾，细胞特异的基因排到前面。

```r
# Geneformer 的输入：把每个细胞变成“按相对表达量排序的基因序列”
cell   <- c(GAPDH = 120, ACTB = 90, MMP9 = 30, IL1B = 15, ACP5 = 8, CTSK = 5, TET2 = 2, OSM = 0)
ref_md <- c(GAPDH = 100, ACTB = 80, MMP9 = 5, IL1B = 4, ACP5 = 1, CTSK = 1, TET2 = 1, OSM = 1)   # 假想的全语料库中位数
rel    <- cell / ref_md                                 # 除以该基因在全语料库里的典型水平，压低“管家基因”
tokens <- names(sort(rel[rel > 0], decreasing = TRUE))
tokens
tokens_ko <- setdiff(tokens, "ACP5")                    # “计算机里敲除 ACP5”：把它从序列里删掉
rbind(原始 = match(c("MMP9", "CTSK", "TET2"), tokens), 删除ACP5后 = match(c("MMP9", "CTSK", "TET2"), tokens_ko))
```

```text
[1] "ACP5"  "MMP9"  "CTSK"  "IL1B"  "TET2"  "GAPDH" "ACTB" 
           [,1] [,2] [,3]
原始          2    3    5
删除ACP5后    1    2    4
```

序列里 ACP5 排第一，管家基因排在最后。"计算机里敲除 ACP5"就是把它从序列里删掉：其余基因的顺位整体上移（MMP9 从第 2 变第 1，CTSK 从第 3 变第 2，TET2 从第 5 变第 4）。Geneformer 再把删前和删后的细胞送进模型，看细胞的表征向"疾病"或"健康"状态移动了多少。

### 9.6 怎么回测，怎么对照

报告说要"用已知答案做回测"。回测的做法：

1. **选已有定论的扰动。** 如 TET2 敲除、MMP9 抑制、OSM 中和、IL-20 中和；
2. **定义指标。** 常用的有三个：
   - 方向一致率：预测的评分变化方向与已发表结果是否一致；
   - 排序能力：已知的下游基因在预测排序里的 AUC（见上面的玩具）；
   - 与线性基线的差值：复杂模型比基线好多少；
3. **并排报告基线。** Ahlmann-Eltze 等（Nat Methods 2025）发现，当时的深度学习扰动预测并没有跑赢简单的线性基线；报告引用的 Nat Methods 2025 一篇也得出了相同的结论。所以任何复杂模型的结果都必须同时给出基线；
4. **随机对照。** 随机选同样数量的基因做"虚拟敲除"，看评分变化的分布，给真实候选节点一个经验 p 值。

### 9.7 "双向判定"的逻辑和陷阱

报告的判定是：候选节点在髓系细胞里被抑制后，应当**同时**降低骨髓里的破骨评分和斑块里的钙化评分。

要注意：

1. **两个组织各用各的网络。** 调控网络是细胞类型特异的，不能把骨髓的网络直接套在斑块上。
2. **避免循环论证。** 如果评分用到的基因里包含候选节点自己（如 MMP9 既是候选节点又是评分基因），敲除它评分必然下降。评分基因集要排除候选节点。
3. **效应要大于噪声。** 对照随机节点的分布，别只看符号。
4. **报告两侧的效应量，不只是方向。** "同时降低"要有量化的标准。

### 小结

- CellOracle 借助底图网络、scTenifoldKnk 靠相关性网络、NicheNet 排序配体、Geneformer 靠预训练模型。
- 线性基线经常已经够好；复杂模型必须并排报告基线。
- 方向正确的先验网络是 CellOracle 类方法成败的关键。
- 用已知答案回测，并用随机节点做对照。

---

[← 上一篇：07 · 基因程序评分与跨组织比对](07-gene-program-scoring.md) · [目录](README.md) · [下一篇：09 · 个体级队列分析：CHIP × 骨密度 × 骨折，以及共病 →](09-ukb-cohort-and-comorbidity.md)
