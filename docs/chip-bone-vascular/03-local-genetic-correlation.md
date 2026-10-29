[← 上一篇：02 · 遗传相关：LDSC 与 HDL](02-genetic-correlation-ldsc-hdl.md) · [目录](README.md) · [下一篇：04 · 共定位：coloc 与 SuSiE-coloc →](04-colocalization-coloc-susie.md)

---

# 03 · 局部遗传相关：LAVA 与 SUPERGNOVA

> 这一章带走什么：为什么全局遗传相关等于 0 不代表"没有共享"；局部方法的思路；一个"区域版"的玩具估计，以及它能看见什么。

### 4.1 全局为零，可能是正负抵消

报告引用的 Bone 2026 说：CAC 与 eBMD **没有全局遗传相关**，但有 211 个共享基因，集中在 17 号染色体。

这并不矛盾。全局 rg 是整个基因组"平均"出来的一个数。设想基因组分成两半：

- 一半里，让性状 1 升高的变异也让性状 2 升高（局部相关为 +0.9）；
- 另一半里，反过来（局部相关为 −0.9）。

两半一平均，全局 rg 约等于 0。但这两个性状明明有强烈的、结构化的共享。如果还只看全局数，就会得出"没有关系"的错误结论。

### 4.2 局部方法的思路

局部遗传相关方法（报告提到 LAVA 和 SUPERGNOVA）的做法是：

1. 把基因组切成约 2,500 个近似独立的 LD 区块（每块约 1 Mb，具体划分以原文为准）；
2. 在每个区块里，分别估计两个性状的局部遗传度，并检验是否显著。**没有局部遗传信号，就谈不上局部相关**；
3. 只在两个性状局部遗传度都显著的区块里，估计局部遗传协方差和局部 rg，给出 p 值；
4. 样本重叠用 LDSC 的跨性状截距校正（见 [02 遗传相关](02-genetic-correlation-ldsc-hdl.md)）；
5. 区块数量大，要做多重检验校正（如 Bonferroni，阈值约 0.05 ÷ 检验的区块数）。

LAVA 和 SUPERGNOVA 做的是同一件事，估计量的细节不同。报告写"二选一"，我的建议是两个都跑，互相当交叉验证。

与全局方法相比，局部方法的代价是：每个区块的 SNP 少、信号弱，所以**统计功效更低、标准误更大**。

### 4.3 实战：区域版估计

我们沿用 [02 遗传相关](02-genetic-correlation-ldsc-hdl.md) 里的迷你基因组。如果你是直接打开本篇、还没有在当前 R 会话里跑过上一篇，先运行下面这段前置代码把迷你基因组和辅助函数建好（已经跑过上一篇的可以跳过这一段）：

```r
set.seed(2026)
M <- 4000; N <- 8000                                 # 4000 个 SNP；每个 GWAS 样本 8000 人
bsizes <- c()                                        # 连锁块大小不一：5、10、20、40、80 个 SNP
while (sum(bsizes) < M) bsizes <- c(bsizes, sample(c(5, 10, 20, 40, 80), 1))
bsizes[length(bsizes)] <- bsizes[length(bsizes)] - (sum(bsizes) - M)
blk    <- rep(seq_along(bsizes), bsizes)             # 每个 SNP 属于哪个连锁块
rho_b  <- runif(length(bsizes), 0.2, 0.9)            # 每个块内部的相关强度
maf    <- runif(M, 0.10, 0.50)                       # 效应等位基因频率
load_s <- rnorm(M, 0, 0.04)                          # 每个 SNP 受“祖源因子”影响的程度（用来制造混杂）

sim_people <- function(n) {
  anc <- rnorm(n)                                    # 看不见的祖源因子
  G <- matrix(0L, n, M)
  for (b in seq_along(bsizes)) {
    cols <- which(blk == b)
    u <- rnorm(n)                                    # 块内共享的成分，决定块内 SNP 一起变
    Z <- sqrt(rho_b[b]) * u + sqrt(1 - rho_b[b]) * matrix(rnorm(n * length(cols)), n)
    Z <- Z * rep(sqrt(1 - load_s[cols]^2), each = n) + outer(anc, load_s[cols])
    G[, cols] <- (Z > rep(qnorm((1 - maf[cols])^2), each = n)) +     # 阈值切成 0/1/2
                 (Z > rep(qnorm(1 - maf[cols]^2), each = n))
  }
  list(G = G, anc = anc)
}

s1 <- sim_people(N); Gs1 <- scale(s1$G); s1$G <- NULL        # GWAS 样本 1（标准化基因型）
s2 <- sim_people(N); Gs2 <- scale(s2$G); s2$G <- NULL        # GWAS 样本 2，与样本 1 是不同的人
ref <- sim_people(3000)                                      # 另一批人，专门用来算 LD 分数

# LD 分数：每个 SNP 与块内所有 SNP（含自己）的校正 r² 之和
ld <- numeric(M)
for (b in seq_along(bsizes)) {
  cols <- which(blk == b)
  r2 <- cor(ref$G[, cols, drop = FALSE])^2
  ld[cols] <- rowSums(r2 - (1 - r2) / (3000 - 2))
}
rm(ref)
c(SNP数 = M, 连锁块数 = length(bsizes))
round(summary(ld), 1)

gwas_z <- function(Gs, y) {                          # 迷你 GWAS：每个 SNP 一次回归，返回 z 分数
  n <- length(y); ys <- as.vector(scale(y))
  r <- as.vector(crossprod(Gs, ys)) / (n - 1)
  r * sqrt((n - 2) / (1 - r^2))
}

ldsc_uni <- function(z, N, idx) {                    # 单性状：chi2 对 LD 分数做加权回归
  y <- z[idx]^2; l <- ld[idx]
  h0 <- max(unname(coef(lm(y ~ l))[2]) * M / N, 0.001)
  f <- lm(y ~ l, weights = 1 / (l * (1 + N * h0 * l / M)^2))
  c(int = unname(coef(f)[1]), h2 = unname(coef(f)[2]) * M / N)
}

ldsc_cross <- function(z1, z2, N1, N2, idx) {        # 两性状：z1*z2 对 LD 分数做加权回归
  u1 <- ldsc_uni(z1, N1, idx); u2 <- ldsc_uni(z2, N2, idx)
  y <- z1[idx] * z2[idx]; l <- ld[idx]
  w <- 1 / (l * (1 + N1 * u1[["h2"]] * l / M) * (1 + N2 * u2[["h2"]] * l / M))
  f <- lm(y ~ l, weights = w)
  rho_g <- unname(coef(f)[2]) * M / sqrt(N1 * N2)    # 遗传协方差
  c(h2_1 = u1[["h2"]], h2_2 = u2[["h2"]], int_1 = u1[["int"]], int_2 = u2[["int"]],
    int_x = unname(coef(f)[1]), rho_g = rho_g, rg = rho_g / sqrt(max(u1[["h2"]] * u2[["h2"]], 1e-6)))
}

jackknife <- function(fun, nj = 40) {                # 块刀切：每次去掉 1/40 的基因组重算
  full <- fun(seq_len(M)); grp <- cut(seq_len(M), nj, labels = FALSE)
  jk <- t(sapply(seq_len(nj), function(i) fun(which(grp != i))))
  list(est = full, se = sqrt((nj - 1) / nj * colSums(sweep(jk, 2, colMeans(jk))^2)))
}
```

有了迷你基因组，设计如下：

- 把基因组分成 6 个区域（各含约 1/6 的 SNP）；
- 前 3 个区域里，两个性状的效应相关 +0.9；后 3 个区域里相关 −0.9；
- 两个性状来自两批不重叠的人。

下面是一个"区域版"的估计方法，用的还是 LDSC 的思路：在每个区域内，用 LD 分数加权，估计局部遗传度和局部协方差，再相除。**这是教学简化版，不是 LAVA 本身**（LAVA 用的是区块内完整的 LD 矩阵）。为了让演示干净，截距直接取真值（没有混杂、没有重叠：1、1、0）。

```r
set.seed(21)
mid <- cumsum(bsizes) - bsizes / 2                    # 每个连锁块在基因组上的中点位置
region_of_block <- cut(mid, 6, labels = FALSE)        # 把基因组分成 6 个区域，各含约 1/6 的 SNP
region <- region_of_block[blk]                        # 每个 SNP 属于哪个区域
rho_r  <- c(0.9, 0.9, 0.9, -0.9, -0.9, -0.9)          # 前 3 个区域效应同向，后 3 个区域效应反向
h2_loc <- 0.5

two_traits <- function(rho_vec) {                     # 造两个性状（样本不重叠）的 z 分数
  b1 <- rnorm(M, 0, sqrt(h2_loc / M))
  b2 <- rho_vec[region] * b1 + sqrt(1 - rho_vec[region]^2) * rnorm(M, 0, sqrt(h2_loc / M))
  list(z1 = gwas_z(Gs1, as.vector(Gs1 %*% b1) + rnorm(N, 0, sqrt(1 - h2_loc))),
       z2 = gwas_z(Gs2, as.vector(Gs2 %*% b2) + rnorm(N, 0, sqrt(1 - h2_loc))))
}

# 区域版估计（教学简化版，不是 LAVA 本身）：
#   局部协方差 = sum(z1*z2*LD分数) / (N*sum(LD分数^2))；局部遗传度 = sum((z^2-1)*LD分数) / (N*sum(LD分数^2))
#   局部 rg = 局部协方差 / sqrt(局部遗传度1 * 局部遗传度2)
# 截距直接取真值（没有混杂、没有重叠：1、1、0）
local_rg <- function(z1, z2, blocks_in) {
  idx <- which(blk %in% blocks_in); l <- ld[idx]; d <- N * sum(l^2)
  v1 <- sum((z1[idx]^2 - 1) * l) / d
  v2 <- sum((z2[idx]^2 - 1) * l) / d
  sum(z1[idx] * z2[idx] * l) / d / sqrt(max(v1 * v2, 1e-12))
}
```

先看每个区域的构成：

```r
data.frame(区域 = 1:6, SNP数 = as.vector(table(region)), 连锁块数 = as.vector(table(region_of_block)),
           平均LD分数 = round(tapply(ld, region, mean), 1))
```

```text
  区域 SNP数 连锁块数 平均LD分数
1    1   665       19       14.6
2    2   730       16       15.6
3    3   585       17        7.8
4    4   715       20       12.4
5    5   605       15        7.2
6    6   700       23       17.2
```

各区域含 585 到 730 个 SNP、15 到 23 个连锁块，平均 LD 分数在 7 到 17 之间。

**重复 20 次**，比较全局方法和区域方法。另外加一个对照：6 个区域都同向。

```r
nrep <- 20
loc <- matrix(NA, nrep, 6); glob <- numeric(nrep); glob_same <- numeric(nrep)
for (k in 1:nrep) {
  w <- two_traits(rho_r)                                                        # 正负区域各半
  glob[k] <- ldsc_cross(w$z1, w$z2, N, N, seq_len(M))[["rg"]]
  loc[k, ] <- sapply(1:6, function(r) local_rg(w$z1, w$z2, which(region_of_block == r)))
  w2 <- two_traits(rep(0.9, 6))                                                 # 对照：6 个区域都同向
  glob_same[k] <- ldsc_cross(w2$z1, w2$z2, N, N, seq_len(M))[["rg"]]
}
cat("全局 LDSC rg（重复", nrep, "次）：正负各半 平均", round(mean(glob), 2), "SD", round(sd(glob), 2),
    "；全部同向 平均", round(mean(glob_same), 2), "SD", round(sd(glob_same), 2), "\n")
out <- rbind(真实相关 = rho_r, 局部rg平均 = colMeans(loc), 局部rgSD = apply(loc, 2, sd))
colnames(out) <- paste0("区域", 1:6)
round(out, 2)
```

```text
全局 LDSC rg（重复 20 次）：正负各半 平均 -0.09 SD 0.28 ；全部同向 平均 0.91 SD 0.05 
           区域1 区域2 区域3 区域4 区域5 区域6
真实相关    0.90  0.90  0.90 -0.90 -0.90 -0.90
局部rg平均  0.83  0.88  0.92 -0.89 -0.91 -0.89
局部rgSD    0.24  0.09  0.10  0.22  0.08  0.18
```

读法：

- **全局 LDSC：** 正负各半时，rg 平均 −0.09，SD 0.28，完全看不出关系；而 6 个区域都同向时，rg 平均 0.91，SD 0.05。说明全局估计器本身没问题，是正负抵消把它蒙住了。
- **区域版：** 6 个区域的局部 rg 平均分别是 0.83、0.88、0.92、−0.89、−0.91、−0.89，和真值 ±0.9 吻合。符号和大小都对。

**实际分析只有一份数据**，用逐块刀切给出每个区域的标准误：

```r
# 实际分析只有一份数据：逐块刀切给出每个区域的标准误
w <- two_traits(rho_r)
res <- t(sapply(1:6, function(r) {
  bl  <- which(region_of_block == r)
  est <- local_rg(w$z1, w$z2, bl)
  jk  <- sapply(seq_along(bl), function(i) local_rg(w$z1, w$z2, bl[-i]))
  se  <- sqrt((length(bl) - 1) / length(bl) * sum((jk - mean(jk))^2))
  c(真实相关 = rho_r[r], 局部rg = est, SE = se, z = est / se)
}))
rownames(res) <- paste0("区域", 1:6)
round(res, 2)
```

```text
      真实相关 局部rg   SE      z
区域1      0.9   0.92 0.08  11.20
区域2      0.9   0.90 0.21   4.38
区域3      0.9   0.67 0.24   2.77
区域4     -0.9  -0.88 0.09  -9.91
区域5     -0.9  -0.96 0.03 -31.91
区域6     -0.9  -0.90 0.09  -9.57
```

6 个区域的符号全部正确。区域 3 的估计 0.67（SE 0.24）偏离真值 0.9 约 1 个 SE，区域 2 的 SE 也较大（0.21），这是局部估计"信号弱、噪声大"的典型表现。

### 4.4 使用时的要点

1. **先过滤再检验。** 只在两个性状局部遗传度都显著的区块里估计局部 rg，否则 rg 的分母接近 0，结果不稳。CHIP 的样本量小，很多区块会被过滤掉，这是功效问题，不是方法问题。
2. **多重检验。** 区块数约 2,500，每对性状都要校正；报告里 4 个 CHIP 表型 × 3 个结局表型共 12 对，校正的强度要写在方案里（[11 流水线与立项规则](11-pipeline-and-go-nogo.md) 会算这条的后果）。
3. **局部相关不是因果。** 它只说明"这一段基因组上两个性状的遗传效应共变"。接下来要用共定位判断是不是同一个变异（[04 共定位](04-colocalization-coloc-susie.md)）。
4. **方向要两头看。** 正相关和负相关的区块都要报告，不要只挑一种。

### 小结

- 全局 rg 是平均数，正负抵消会把它压成 0。
- 局部方法逐区块估计，只在有局部遗传信号的区块里检验。
- 代价是功效更低、要做多重检验；CHIP 的有效样本量小，需要有心理准备。

---

[← 上一篇：02 · 遗传相关：LDSC 与 HDL](02-genetic-correlation-ldsc-hdl.md) · [目录](README.md) · [下一篇：04 · 共定位：coloc 与 SuSiE-coloc →](04-colocalization-coloc-susie.md)
