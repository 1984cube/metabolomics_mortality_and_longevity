[← 上一篇：04 · 共定位：coloc 与 SuSiE-coloc](04-colocalization-coloc-susie.md) · [目录](README.md) · [下一篇：06 · 把 GWAS 信号映射到细胞：MAGMA、scDRS 与 scPagwas →](06-gwas-to-cells-magma-scdrs.md)

---

# 05 · 单细胞分析基础，以及“别把细胞当重复”

> 这一章带走什么：单细胞数据的标准分析流程；为什么单细胞研究的"样本"是供体而不是细胞；把细胞当重复会造成多严重的假阳性；17 个携带者 vs 7 个对照的数据能检出多大的差异。

### 6.1 数据长什么样

单细胞 RNA-seq 的原始结果是一张**计数矩阵**：行是细胞，列是基因，格子里是这个细胞里读到的这个基因的分子数。它很稀疏（大多数格子是 0），而且每个细胞的总读数（测序深度）差别很大。

### 6.2 标准流程

| 步骤 | 做什么 | 常用工具 |
|---|---|---|
| 质控 | 去掉低质量细胞、双细胞、线粒体比例过高的细胞 | Seurat（R）、Scanpy（Python） |
| 标准化 | 按每个细胞的总读数校正深度，再取 log1p | 同上 |
| 高变基因 | 选出变化最大的一批基因，降低噪声 | 同上 |
| 降维 | 主成分分析（PCA）；可视化用 UMAP | 同上 |
| 聚类 | 把相似的细胞分成群 | 同上 |
| 注释 | 用标记基因给每个群命名（T 细胞、巨噬细胞……） | 同上 |
| 差异分析 | 比较组间的基因表达；**在供体层面做** | edgeR、DESeq2、limma-voom |

### 6.3 实战 1：用基础 R 走一遍流程

造一个玩具数据：3 种细胞类型各 500 个，600 个基因，每种类型有 50 个标记基因表达更高，每个细胞的测序深度不同。

```r
set.seed(31)
ngene <- 600; per <- 500
type  <- rep(c("T细胞", "B细胞", "髓系"), each = per); ncell <- length(type)
base  <- rnorm(ngene, 0, 1.0)                             # 每个基因的基础表达量（对数尺度）
fc    <- matrix(0, ngene, 3)                              # 每种细胞类型各有 50 个标记基因表达更高
for (k in 1:3) fc[sample(ngene, 50), k] <- rnorm(50, 1.1, 0.3)
depth <- rlnorm(ncell, 0, 0.4)                            # 每个细胞的测序深度不同
lambda <- depth * exp(matrix(base, ncell, ngene, byrow = TRUE) + t(fc)[match(type, unique(type)), ])
counts <- matrix(rpois(length(lambda), lambda), ncell)    # 计数矩阵：行是细胞，列是基因

lib <- rowSums(counts)
X   <- log1p(counts / lib * 1e4)                          # 1. 按测序深度标准化，再取 log1p
hvg <- order(apply(X, 2, var), decreasing = TRUE)[1:200]  # 2. 选 200 个变化最大的基因
pcs <- prcomp(X[, hvg], scale. = TRUE)$x[, 1:10]          # 3. 主成分分析，取前 10 个
cl  <- kmeans(pcs, 3, nstart = 10)$cluster                # 4. 聚类
table(真实类型 = type, 聚类 = cl)
```

```text
        聚类
真实类型   1   2   3
   B细胞 500   0   0
   T细胞   0 500   0
   髓系    0   0 500
```

聚类结果和真实类型完全一致。**这个玩具数据特意做得很干净**，真实数据里亚型之间界线模糊，聚类不会这么整齐。

再找每个聚类的标记基因，看能不能找回真实的标记基因：

```r
# 5. 找每个聚类的标记基因：该聚类里的均值减去其余细胞的均值，取前 20 个
true_marker <- lapply(1:3, function(k) which(fc[, k] > 0))
res <- lapply(1:3, function(k) {
  top  <- order(colMeans(X[cl == k, ]) - colMeans(X[cl != k, ]), decreasing = TRUE)[1:20]
  best <- which.max(sapply(1:3, function(t) mean(top %in% true_marker[[t]])))   # 这个聚类最像哪种真实类型
  data.frame(聚类 = k, 对应类型 = unique(type)[best], 前20个里真标记基因的个数 = sum(top %in% true_marker[[best]]))
})
do.call(rbind, res)
```

```text
  聚类 对应类型 前20个里真标记基因的个数
1    1    B细胞                       20
2    2    T细胞                       20
3    3     髓系                       20
```

三个聚类各自对应一种真实类型，各自的前 20 个标记基因全部是真标记基因。

### 6.4 报告里的数据集分级

报告把单细胞数据集分成 A、B、C 三级：

- **A：可直接做组间比较。** 看的是**每组的供体数**，不是细胞数。例如 CNGB TrueBlood 有 17 例携带者 + 7 例对照；GSE159677 虽然只有 3 对配对样本，但配对设计很关键，也归 A；
- **B：可做状态发现与方向验证。** 样本少，不适合严格的组间检验；
- **C：参考用。**

### 6.5 最大的坑：伪重复

设想 CNGB 数据：17 例携带者 + 7 例对照，每人取 200 个细胞，共 4,800 个细胞。你想找"携带者 vs 对照"表达有差别的基因。

直觉做法：把 4,800 个细胞分成两组，每个基因做一次 t 检验。**这是错的。** 同一个供体的 200 个细胞共享这个供体自己的"基线"（遗传背景、年龄、批次、采样当天的状态……），它们不是 200 个独立样本，更像是对同一个人的 200 次重复测量。

模拟验证：1,000 个基因，其中只有前 50 个基因真的有组间差异，其余 950 个没有。供体之间有天然差异（SD 0.35）。

```r
set.seed(32)
n_case <- 17; n_ctrl <- 7; cells <- 200; ngene <- 1000      # 17 例携带者 + 7 例对照，每人取 200 个细胞
grp    <- rep(c("携带者", "对照"), c(n_case, n_ctrl))
donor  <- rep(seq_along(grp), each = cells); g_cell <- grp[donor]
base   <- rnorm(ngene, 0.3, 1)
true_fc <- c(rep(0.5, 50), rep(0, ngene - 50))              # 前 50 个基因有真实的组间差异；其余 950 个没有
donor_eff <- matrix(rnorm(length(grp) * ngene, 0, 0.35), length(grp), ngene)   # 供体之间的天然差异
logmu  <- matrix(base, length(donor), ngene, byrow = TRUE) + donor_eff[donor, ] +
          outer(as.integer(g_cell == "携带者"), true_fc)
depth  <- rlnorm(length(donor), 0, 0.3)
counts <- matrix(rpois(length(logmu), depth * exp(logmu)), length(donor))
X <- log1p(counts / rowSums(counts) * mean(rowSums(counts)))

welch_p <- function(M, is1) {                                # 对每一列做 Welch t 检验
  m1 <- colMeans(M[is1, , drop = FALSE]);  m0 <- colMeans(M[!is1, , drop = FALSE])
  v1 <- apply(M[is1, , drop = FALSE], 2, var) / sum(is1); v0 <- apply(M[!is1, , drop = FALSE], 2, var) / sum(!is1)
  t <- (m1 - m0) / sqrt(v1 + v0); df <- (v1 + v0)^2 / (v1^2 / (sum(is1) - 1) + v0^2 / (sum(!is1) - 1))
  2 * pt(-abs(t), df)
}
p_cell <- welch_p(X, g_cell == "携带者")                      # 做法 A：把每个细胞当一个独立样本
pb     <- rowsum(X, donor) / cells                           # 做法 B：先按供体求平均（伪批量），再比较 24 个供体
p_pb   <- welch_p(pb, grp == "携带者")

summ <- function(p) {
  disc <- p.adjust(p, "BH") < 0.05
  c(`无差异基因中 p<0.05 的比例` = mean(p[51:ngene] < 0.05), `真有差异的基因中 p<0.05 的比例` = mean(p[1:50] < 0.05),
    `FDR<5% 的发现数` = sum(disc), `其中假发现` = sum(disc[51:ngene]))
}
round(rbind(`细胞当重复` = summ(p_cell), `按供体汇总` = summ(p_pb)), 2)
```

```text
           无差异基因中 p<0.05 的比例 真有差异的基因中 p<0.05 的比例 FDR<5% 的发现数 其中假发现
细胞当重复                       0.67                           1.00             672        623
按供体汇总                       0.05                           0.88               8          1
```

| 做法 | 无差异的 950 个基因里，p<0.05 的比例 | FDR<5% 的发现数 | 其中假发现 |
|---|---|---|---|
| 细胞当重复 | 67% | 672 | 623（93%） |
| 按供体汇总 | 5% | 8 | 1 |

**把细胞当重复，三分之二的无关基因都"显著"，报出来的发现里 93% 是假的。** 按供体汇总后，假阳性率回到 5%，发现很少但可信。细胞当重复的"灵敏度 100%"是假象：它对一切都敏感。

### 6.6 怎么做对

三种常用做法：

1. **伪批量（pseudobulk）。** 把同一个供体、同一种细胞类型的原始计数加起来，得到每个供体一份"表达谱"，再用 edgeR、DESeq2 或 limma-voom 这类工具比较两组供体。可以方便地加入年龄、性别、批次等协变量。
2. **混合模型。** 保留细胞层面的数据，把供体设为随机效应（如 lme4、glmmTMB、NEBULA）。信息利用更充分，计算更重。
3. **供体层面的置换检验。** 打乱供体的组别标签，重新算统计量。

三种做法的共同点：**重复单位是供体**。报告结果时要写明每组的供体数。参考：Squair 等，Nat Commun 2021；Zimmerman 等，Nat Commun 2021。

### 6.7 17 vs 7 个供体，能检出多大的差异

样本量小，要先知道"功效够不够"。下面算 80% 功效、双侧检验时，最小可检出的标准化效应（单位：供体间的标准差）。

```r
# 17 例携带者 vs 7 例对照：80% 功效能检出多大的差异？（单位：供体间标准差）
power_unbal <- function(d, n1, n0, alpha = 0.05) {
  df <- n1 + n0 - 2; ncp <- d * sqrt(n1 * n0 / (n1 + n0)); tc <- qt(1 - alpha / 2, df)
  1 - pt(tc, df, ncp) + pt(-tc, df, ncp)
}
min_d <- function(n1, n0, alpha = 0.05) uniroot(function(d) power_unbal(d, n1, n0, alpha) - 0.8, c(0.1, 6))$root
data.frame(设计 = c("17 vs 7", "17 vs 17", "7 vs 7", "17 vs 7（每个基因 α=0.001）"),
           最小可检出效应 = round(c(min_d(17, 7), min_d(17, 17), min_d(7, 7), min_d(17, 7, 0.001)), 2))
```

```text
                         设计 最小可检出效应
1                     17 vs 7           1.32
2                    17 vs 17           0.99
3                      7 vs 7           1.63
4 17 vs 7（每个基因 α=0.001）           2.12
```

17 vs 7 的设计，最小可检出效应是 1.32 个供体间标准差。也就是说，组间差异要大于一个多标准差，才有 80% 的把握检出。如果要对每个基因做严格的多重检验（α = 0.001），门槛涨到 2.12。**所以这类数据适合发现"大效应"和"程序层面的整体方向"，不适合逐个基因地找小效应。** 报告里把 A 级数据用于"组间比较"、其余用于"方向验证"，正是这个道理。

### 小结

- 单细胞流程：标准化、高变基因、PCA、聚类、注释、差异分析。
- 重复单位是供体，不是细胞；把细胞当重复，假阳性率可以高达三分之二。
- 17 vs 7 的设计只能检出大效应（约 1.3 个供体间标准差以上）。

---

[← 上一篇：04 · 共定位：coloc 与 SuSiE-coloc](04-colocalization-coloc-susie.md) · [目录](README.md) · [下一篇：06 · 把 GWAS 信号映射到细胞：MAGMA、scDRS 与 scPagwas →](06-gwas-to-cells-magma-scdrs.md)
