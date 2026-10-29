[← 上一篇：06 · 把 GWAS 信号映射到细胞：MAGMA、scDRS 与 scPagwas](06-gwas-to-cells-magma-scdrs.md) · [目录](README.md) · [下一篇：08 · 虚拟扰动：CellOracle、scTenifoldKnk、NicheNet 与 Geneformer →](08-in-silico-perturbation.md)

---

# 07 · 基因程序评分与跨组织比对

> 这一章带走什么：什么是"基因程序"和程序评分；基于排名的评分为什么更稳健；评分必须在供体层面比较；怎么判断两个组织里的程序是否重叠，以及一个会让你"假阳性"的共同基线陷阱。

### 8.1 报告里的程序

报告要在 CHIP 外周血单细胞数据里构建"破骨样巨噬细胞程序"评分，用 ACP5、CTSK、MMP9、NFATC1、OSM 这五个核心基因。它们各自的角色：

- **ACP5**：抗酒石酸酸性磷酸酶（TRAP），破骨细胞的经典标记；
- **CTSK**：组织蛋白酶 K，破骨细胞降解骨基质的主要蛋白酶；
- **MMP9**：基质金属蛋白酶 9，降解细胞外基质；
- **NFATC1**：破骨细胞分化的关键转录因子；
- **OSM**：oncostatin M，一种细胞因子（报告说 TET2 突变单核细胞通过它促进瓣膜间质细胞钙化）。

然后在钙化斑块和骨髓两套数据里分别评分，**比较程序重叠度**。

### 8.2 程序评分的几种做法

| 做法 | 思路 | 特点 |
|---|---|---|
| 平均表达减对照 | 程序基因的平均表达，减去表达量匹配的随机基因的平均（Seurat 的 AddModuleScore，思路来自 Tirosh 2016） | 简单；对深度、表达量水平敏感 |
| 基于排名（UCell） | 看程序基因在每个细胞里的表达排名（Andreatta 与 Carmona 2021） | 每个细胞独立计算，不受其他细胞影响；对批次和深度更稳健 |
| 基于排名的曲线下面积（AUCell） | 程序基因是否集中在细胞表达最高的那批基因里（Aibar 2017） | 同样基于排名 |

UCell 的公式：

```text
U = 1 − ( Σ r_i − n(n+1)/2 ) ÷ ( n × maxRank )
```

其中 n 是程序基因个数，r_i 是第 i 个程序基因在这个细胞里的表达排名（1 表示最高；排名超过 maxRank 的按 maxRank + 1 算）。n 个程序基因都排在最前面时 U = 1；都排在 maxRank 之外时 U ≈ 0。

### 8.3 实战：两个组织、两种情形

造玩具数据：每个组织 16 个供体（携带者和对照各 8 个），每个供体 150 个巨噬细胞，1,500 个基因。其中前 40 个基因是"破骨样程序"（含五个核心基因）。设计：

- 每个细胞有一定概率处于"破骨样"状态（携带者平均 30%，对照平均 12%，供体之间还有差异）；
- 破骨样细胞里，程序基因的表达升高；
- **骨髓**：程序 = 基因 1 到 40；
- **斑块，情形 S（共享）**：程序与骨髓相同；
- **斑块，情形 T（不同）**：程序是另外 40 个基因，不含五个核心基因。

```r
set.seed(51)
ngene <- 1500; nd <- 8; cells_per <- 150               # 每组 8 个供体，每个供体 150 个巨噬细胞
genes <- c("ACP5", "CTSK", "MMP9", "NFATC1", "OSM", sprintf("PG%02d", 1:35), sprintf("G%04d", 1:(ngene - 40)))
base  <- rnorm(ngene, 0.4, 1.2)                        # 基因的基础表达量（对数尺度）

sim_tissue <- function(program) {                      # program：这个组织里“破骨样程序”包含哪些基因（下标）
  group  <- rep(c("携带者", "对照"), each = nd)
  p_ocl  <- ifelse(group == "携带者", 0.30, 0.12) * exp(rnorm(2 * nd, 0, 0.2))   # 每个供体里破骨样细胞的比例
  donor  <- rep(seq_len(2 * nd), each = cells_per)
  ocl    <- rbinom(length(donor), 1, p_ocl[donor])     # 每个细胞是否处于破骨样状态
  d_eff  <- matrix(rnorm(2 * nd * ngene, 0, 0.15), 2 * nd, ngene)                 # 供体间的天然差异
  logmu  <- matrix(base, length(donor), ngene, byrow = TRUE) + d_eff[donor, ]
  logmu[, program] <- logmu[, program] + 1.5 * ocl     # 破骨样细胞里，程序基因表达升高
  counts <- matrix(rpois(length(logmu), exp(logmu)), length(donor))
  X <- log1p(counts / rowSums(counts) * 1e4); colnames(X) <- genes
  list(X = X, donor = donor, group = group, ocl = ocl)
}
bm     <- sim_tissue(1:40)                             # 骨髓：程序 = 基因 1–40（含 5 个核心基因）
pl_sh  <- sim_tissue(1:40)                             # 斑块（情形 S）：与骨髓共享同一套程序
pl_sp  <- sim_tissue(101:140)                          # 斑块（情形 T）：另一套完全不同的程序
```

**评分和比较。** 用一个 UCell 式的排名评分（教学简化版，细节如并列排名的处理可能与 UCell 软件包不同）给每个细胞打分，再在**供体层面**比较携带者和对照；同时也报告"把细胞当重复"的 p 值，看差别。

```r
# UCell 式的基于排名的评分（教学简化版）：核心基因在每个细胞里的表达排名越靠前，得分越高
ucell <- function(X, sig, maxRank = 300) {
  idx <- match(sig, colnames(X)); n <- length(idx)
  apply(X, 1, function(x) {
    r <- pmin(rank(-x, ties.method = "average")[idx], maxRank + 1)
    1 - (sum(r) - n * (n + 1) / 2) / (n * maxRank)
  })
}
core5 <- c("ACP5", "CTSK", "MMP9", "NFATC1", "OSM")
welch <- function(v, g) t.test(v[g == "携带者"], v[g == "对照"])$p.value

report <- function(tis) {
  sc <- ucell(tis$X, core5)
  d_score <- tapply(sc, tis$donor, mean)                       # 先按供体汇总
  c(`破骨样细胞比例` = mean(tis$ocl),
    `细胞当重复的p` = welch(sc, tis$group[tis$donor]),
    `按供体汇总的p` = welch(d_score, tis$group))
}
signif(rbind(骨髓 = report(bm), 斑块_共享程序 = report(pl_sh), 斑块_不同程序 = report(pl_sp)), 2)
```

```text
              破骨样细胞比例 细胞当重复的p 按供体汇总的p
骨髓                    0.19       2.7e-11       0.01100
斑块_共享程序           0.22       1.4e-28       0.00011
斑块_不同程序           0.21       1.1e-04       0.21000
```

读法：

- **骨髓和斑块（共享程序）：** 供体层面的 p 值分别为 0.011 和 0.00011，显著；
- **斑块（不同程序）：** 核心五基因不在这个组织的程序里，供体层面 p = 0.21，不显著，符合设计；
- 但把细胞当重复时，斑块（不同程序）的 p 值是 1.1e-4，**假阳性**。这又是 [05 单细胞基础](05-single-cell-basics-pseudobulk.md) 讲的伪重复。

### 8.4 跨组织比对：两个组织的程序是否重叠

做法：

1. 在每个组织里，定义"程序"= 破骨样细胞 vs 其余巨噬细胞的差异表达（每个基因的 log 倍数变化）。真实数据里要靠聚类和注释确定哪些是破骨样细胞，这里用真实标签；
2. 比较两个组织的程序：
   - 所有基因的 log 倍数变化的 Spearman 相关；
   - 各取上调最多的前 50 个基因，看重叠，用超几何检验；
3. 用置换检验给相关性一个 p 值。

**陷阱：** 两个组织共用同一批基因，基因的基线表达量有高有低。标准化（按总读数缩放）会让高表达基因的 log 倍数系统性地偏移，这个偏移在两个组织里是**共同的**，会造出一个假的相关。如果置换检验只是随机打乱基因，就看不出来。正确做法是**按表达量分箱后，只在同一档内打乱基因**，保留基线表达量带来的共同结构。

```r
# 程序比对：每个组织里，破骨样细胞 vs 其余细胞的 log 倍数变化（用真实标签；真实数据里要靠聚类和注释）
logfc <- function(tis) colMeans(tis$X[tis$ocl == 1, ]) - colMeans(tis$X[tis$ocl == 0, ])
fc_bm <- logfc(bm)
bin   <- cut(colMeans(bm$X), quantile(colMeans(bm$X), 0:20 / 20), include.lowest = TRUE, labels = FALSE)  # 按平均表达量分 20 档
shuffle_within_bins <- function(v) { out <- v; for (b in unique(bin)) { i <- which(bin == b); out[i] <- v[i][sample(length(i))] }; out }

compare <- function(fc2, nperm = 2000) {
  rho  <- cor(fc_bm, fc2, method = "spearman")
  top1 <- order(fc_bm, decreasing = TRUE)[1:50]; top2 <- order(fc2, decreasing = TRUE)[1:50]
  ov   <- length(intersect(top1, top2))                                      # 各取前 50 个上调基因，看重叠
  p_ov <- phyper(ov - 1, 50, ngene - 50, 50, lower.tail = FALSE)             # 超几何检验
  p1 <- (1 + sum(replicate(nperm, cor(fc_bm, sample(fc2), method = "spearman")) >= rho)) / (nperm + 1)           # 简单置换
  p2 <- (1 + sum(replicate(nperm, cor(fc_bm, shuffle_within_bins(fc2), method = "spearman")) >= rho)) / (nperm + 1)  # 分箱置换
  c(Spearman相关 = rho, 前50重叠 = ov, 重叠p = p_ov, 简单置换p = p1, 分箱置换p = p2)
}
signif(rbind(共享程序 = compare(logfc(pl_sh)), 不同程序 = compare(logfc(pl_sp))), 2)
```

```text
         Spearman相关 前50重叠   重叠p 简单置换p 分箱置换p
共享程序         0.35       40 1.2e-59     5e-04     0.001
不同程序         0.24        0 1.0e+00     5e-04     0.870
```

读法：

| 情形 | Spearman | 前 50 重叠 | 超几何 p | 简单置换 p | 分箱置换 p |
|---|---|---|---|---|---|
| 共享程序 | 0.35 | 40 | 1.2e-59 | 5e-04 | 0.001 |
| 不同程序 | 0.24 | 0 | 1.0 | 5e-04 | 0.87 |

- 共享程序：重叠 40 个基因，三种检验都显著；
- 不同程序：前 50 个基因**一个也没有重叠**，超几何检验不显著，分箱置换也不显著（p = 0.87）。但 Spearman 相关有 0.24，**简单置换 p 却是 5e-04（最小可能值）**，被共同的基线陷阱骗了。

所以报告里"比较程序重叠度"，**不要只报 Spearman 相关，要同时报前 N 个基因的重叠和分箱置换的 p 值**。

### 8.5 其他要注意的点

1. **不要直接比较不同数据集的原始分数。** 平台、深度、批次都不同。比较"方向"（携带者是否高于对照）和"排序"。
2. **程序基因在两个数据集里都要被检出。** 如果某个核心基因在斑块数据里几乎测不到，它就不应该进入评分。
3. **细胞注释的质量决定一切。** "破骨样细胞"是怎么在两个组织里定义的，要写清楚，最好用同一套标准。
4. **假说不成立也是结果。** 如果两个组织的程序不重叠（情形 T），这回答了问题：程序是组织特异的，不是共享的。报告提醒要把论文标题提前写成中性。

### 小结

- 程序评分用基于排名的方法更稳健；在供体层面比较，不要把细胞当重复。
- 跨组织比对要同时报相关、前 N 重叠，以及分箱置换的 p 值。
- 共同的基线表达会造出假的跨组织相关；简单置换看不出来。

---

[← 上一篇：06 · 把 GWAS 信号映射到细胞：MAGMA、scDRS 与 scPagwas](06-gwas-to-cells-magma-scdrs.md) · [目录](README.md) · [下一篇：08 · 虚拟扰动：CellOracle、scTenifoldKnk、NicheNet 与 Geneformer →](08-in-silico-perturbation.md)
