[← 上一篇：01 · CHIP 与 VAF：数据从哪来](01-chip-and-vaf.md) · [目录](README.md) · [下一篇：03 · 局部遗传相关：LAVA 与 SUPERGNOVA →](03-local-genetic-correlation.md)

---

# 02 · 遗传相关：LDSC 与 HDL

> 这一章带走什么：遗传度、遗传协方差和遗传相关是什么；LD 分数回归怎么靠"信号随 LD 分数上升、混杂不随 LD 分数上升"把两者分开；样本重叠为什么只改截距；CHIP 的有效样本量为什么是个大问题。

### 3.1 先认识三个概念

- **遗传度 h²：** 一个性状的差异里，有多大比例能由全基因组 SNP 的加性效应解释。
- **遗传协方差 ρ_g：** 两个性状的 SNP 效应"一起变化"的程度。如果让某个 SNP 的效应使性状 1 升高，也使性状 2 升高，协方差为正。
- **遗传相关 rg：** 把协方差按两个遗传度标准化：

```text
rg = ρ_g ÷ √(h²₁ × h²₂)
```

可以把 rg 想成两张"信号地图"的重合度：1 表示完全重合，0 表示毫无关系，负数表示方向相反。

### 3.2 汇总统计 + LD：LDSC 的出发点

只有汇总统计，没有个体基因型，怎么估计遗传度？

关键事实：SNP 之间互相连锁（LD）。一个因果变异会让它周围一串 SNP 都带上一点信号。所以：

- 一个 SNP 周围"邻居"越多、关联越强，它就越能"代表"周围的因果信号；
- 这个"代表能力"用 LD 分数衡量：ℓ_j = 该 SNP 与周围所有 SNP 的 r² 之和（包括它自己）。例：块内有 20 个 SNP，彼此 r² 约 0.25，则 ℓ ≈ 1 + 19 × 0.25 = 5.75。

对一个多基因性状（成千上万个小效应变异），**LD 分数大的 SNP，平均的 χ² 就更大**。这条规律是 LDSC 的全部出发点。

### 3.3 LDSC 的核心公式

```text
E[χ²_j] = N h² ℓ_j / M  +  N a  +  1
          └ 斜率 × LD分数 ┘   └ 混杂 ┘
```

- N 是样本量，M 是 SNP 个数，ℓ_j 是第 j 个 SNP 的 LD 分数；
- **斜率** = N h² / M，可以反推出遗传度；
- **截距** = 1 + N a，其中 a 衡量人群分层等混杂。没有混杂时截距是 1。

为什么能把两者分开？因为真正的多基因信号随 LD 分数升高；而混杂（如人群分层）带来的虚假关联，大体上与 LD 分数无关。所以：把 χ² 对 LD 分数做回归，**斜率给信号，截距给混杂**。

对比一个更简单的指标 λGC（χ² 的中位数 ÷ 0.455）：它只能看到总的膨胀，分不清是多基因性还是混杂。

### 3.4 两个性状：遗传相关

把 χ² 换成两个性状 z 分数的乘积：

```text
E[z1_j × z2_j] = √(N1 N2) ρ_g ℓ_j / M  +  ρ N_s / √(N1 N2)
                 └ 斜率：遗传协方差 ┘    └ 截距：样本重叠 ┘
```

- 斜率给遗传协方差，再用两个遗传度标准化得到 rg；
- **样本重叠只改变截距，不改变斜率**。其中 N_s 是两份 GWAS 里重叠的人数，ρ 是两个性状在这些人身上的表型相关；
- 所以 LDSC 估计的 rg 对样本重叠是稳健的，这是它相对 MR 的一个优点（MR 要求样本不重叠，见 [MR 讲义](../mendelian-randomization.md)）。

### 3.5 报告里的第一个问题：有效样本量

CHIP GWAS 是病例对照设计，病例很少。病例对照数据的"有效样本量"是：

```text
有效样本量 = 4 ÷ (1/病例数 + 1/对照数)
```

```r
# 病例对照研究的“有效样本量”：4 / (1/病例 + 1/对照)
eff_n <- function(ncase, nctrl) 4 / (1 / ncase + 1 / nctrl)
cases <- c(10203, 5185, 2041, 4049, 6154)
data.frame(性状 = c("CHIP 总体", "DNMT3A-CH", "TET2-CH", "大克隆", "小克隆"),
           病例 = cases, 对照 = 173918, 有效样本量 = round(eff_n(cases, 173918)))
```

```text
       性状  病例   对照 有效样本量
1 CHIP 总体 10203 173918      38550
2 DNMT3A-CH  5185 173918      20140
3   TET2-CH  2041 173918       8069
4    大克隆  4049 173918      15828
5    小克隆  6154 173918      23775
```

CHIP 总体的有效样本量约 3.9 万；TET2-CH 只有约 8,000。与 BMD GWAS 的几十万人相比，CHIP 这一侧的统计功效很低，所以报告把"遗传层算不出信号"列为第一风险是对的。

经验上，两个性状各自的遗传度估计都要足够准（常用的经验线是遗传度的 z 分数不低于约 4），rg 才可信。具体以 LDSC 文档为准。

### 3.6 实战：搭一个迷你基因组，亲手做 LDSC

真实 LDSC 要用上百万个 SNP 和参考面板。我们造一个迷你版，看原理对不对。

这段代码做了四件事：

1. 造 4000 个 SNP，分成大小不同的连锁块（5 到 80 个 SNP），每块有各自的块内相关强度；
2. 造三批人：两个 GWAS 样本（各 8000 人，互不重叠）和一个只用来算 LD 分数的参考样本（3000 人）；
3. 每个人有一个"祖源因子"，同时轻微影响基因型和表型，用来制造混杂；
4. 算每个 SNP 的 LD 分数。

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
```

```text
   SNP数 连锁块数 
    4000      110 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    1.1     4.0     9.3    12.7    20.2    43.3 
```

迷你基因组只有 110 个连锁块。LD 分数从 1.1 到 43.3，范围够宽，回归才有"杠杆"。

下一步写三个函数：迷你 GWAS（每个 SNP 一次回归）、LDSC 回归、块刀切。

```r
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

**遗传度和混杂的分离。** 造 40 个性状，真实遗传度都是 0.30。一半没有混杂，一半有混杂（祖源因子同时影响基因型和表型）。

```r
set.seed(11)
R <- 40; h2_true <- 0.30                              # 重复造 40 个表型，看估计的平均值和波动
B <- matrix(rnorm(M * R, 0, sqrt(h2_true / M)), M, R) # 每次重复各自的真实效应：每个 SNP 的方差 = h2/M
sim_y <- function(Gs, anc, B, gamma) {                # 表型 = 遗传 + gamma×祖源因子 + 噪声
  n <- nrow(Gs)
  Gs %*% B + gamma * anc + matrix(rnorm(n * ncol(B), 0, sqrt(1 - h2_true - gamma^2)), n)
}
run_h2 <- function(gamma) {
  Y <- sim_y(Gs1, s1$anc, B, gamma)
  t(sapply(seq_len(R), function(i) ldsc_uni(gwas_z(Gs1, Y[, i]), N, seq_len(M))))
}
r0 <- run_h2(0)                                       # 没有混杂
r1 <- run_h2(0.25)                                    # 有混杂（祖源因子同时影响基因型和表型）
data.frame(混杂 = c("无", "有"),
           截距_平均 = round(c(mean(r0[, "int"]), mean(r1[, "int"])), 2),
           截距_SD   = round(c(sd(r0[, "int"]),   sd(r1[, "int"])), 2),
           h2_平均   = round(c(mean(r0[, "h2"]),  mean(r1[, "h2"])), 2),
           h2_SD     = round(c(sd(r0[, "h2"]),    sd(r1[, "h2"])), 2))
```

```text
  混杂 截距_平均 截距_SD h2_平均 h2_SD
1   无       1.0    0.51    0.29  0.06
2   有       1.7    0.49    0.28  0.07
```

读法：

- 没有混杂时，截距平均 1.0，遗传度平均 0.29（真值 0.30）。
- 有混杂时，截距升到平均 1.7，而遗传度估计几乎不变（0.28）。**这就是"截距吸收混杂、斜率只管信号"。**
- 单次估计的波动很大：截距的 SD 约 0.5。原因是迷你基因组只有 110 个连锁块，真实数据有上百万个 SNP，截距的标准误常常只有 0.01 到 0.02。

**实际分析只有一份数据。** 这时用块刀切估计标准误：每次去掉基因组的 1/40，重算一遍，看估计的波动。

```r
# 实际分析时只有一份数据：用块刀切估计标准误
z_one <- gwas_z(Gs1, sim_y(Gs1, s1$anc, B[, 1, drop = FALSE], 0.25)[, 1])
jk <- jackknife(function(idx) ldsc_uni(z_one, N, idx))
round(rbind(估计 = jk$est, 刀切SE = jk$se), 2)
```

```text
        int   h2
估计   2.48 0.20
刀切SE 0.45 0.06
```

这一份数据的截距是 2.48，刀切标准误 0.45，与上面 40 次重复的实际波动（SD 0.49）吻合，说明刀切标准误是可信的。

**遗传相关与样本重叠。** 两个性状的真实 rg = 0.5。重复 30 次，两种设计：

- 不重叠：两个性状来自两批不同的人；
- 完全重叠：两个性状在同一批人身上测，环境噪声相关 0.4。

```r
set.seed(12)
rho_beta <- 0.5                                       # 两个性状的真实效应相关 0.5：真实 rg = 0.5
noise_sd <- sqrt(1 - h2_true)
one_rep <- function() {
  b1 <- rnorm(M, 0, sqrt(h2_true / M))
  b2 <- rho_beta * b1 + sqrt(1 - rho_beta^2) * rnorm(M, 0, sqrt(h2_true / M))
  # 情形一：两个性状来自两批不同的人（样本不重叠）
  zA <- gwas_z(Gs1, as.vector(Gs1 %*% b1) + rnorm(N, 0, noise_sd))
  zB <- gwas_z(Gs2, as.vector(Gs2 %*% b2) + rnorm(N, 0, noise_sd))
  a <- ldsc_cross(zA, zB, N, N, seq_len(M))
  # 情形二：两个性状在同一批人身上测（样本完全重叠），环境噪声相关 0.4
  e  <- matrix(rnorm(2 * N), N, 2) %*% chol(matrix(c(1, 0.4, 0.4, 1), 2)) * noise_sd
  y1 <- as.vector(Gs1 %*% b1) + e[, 1]; y2 <- as.vector(Gs1 %*% b2) + e[, 2]
  o <- ldsc_cross(gwas_z(Gs1, y1), gwas_z(Gs1, y2), N, N, seq_len(M))
  c(rg_不重叠 = a[["rg"]], 截距_不重叠 = a[["int_x"]], rg_重叠 = o[["rg"]], 截距_重叠 = o[["int_x"]],
    表型相关 = cor(y1, y2))
}
reps <- t(replicate(30, one_rep()))
out <- rbind(平均 = colMeans(reps), SD = apply(reps, 2, sd))
round(out, 2)
```

```text
     rg_不重叠 截距_不重叠 rg_重叠 截距_重叠 表型相关
平均      0.48        0.02    0.49      0.39     0.43
SD        0.14        0.26    0.14      0.30     0.01
```

读法：

- 两种设计里 rg 的平均都接近 0.5（0.48 和 0.49）：**样本重叠没有改变 rg**。
- 跨性状截距：不重叠时约 0.02（≈ 0），完全重叠时约 0.39，接近两个性状的表型相关 0.43。**样本重叠只出现在截距里。**

### 3.7 真实使用时的要点

官方 LDSC 是命令行工具（Python）。流程示意如下（依据 LDSC 文档，没有在本环境运行）：

```bash
# 1. 把汇总统计整理成 LDSC 格式（只保留 HapMap3 SNP，统一等位基因）
munge_sumstats.py --sumstats chip.tsv.gz --snp variant_id --a1 effect_allele --a2 other_allele \
    --p p_value --signed-sumstats beta,0 --N-cas 10203 --N-con 173918 \
    --merge-alleles w_hm3.snplist --out chip

# 2. 算两个性状的遗传相关；二分类性状要给出样本中和人群中的患病率，用于换算尺度
ldsc.py --rg chip.sumstats.gz,ebmd.sumstats.gz \
    --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ \
    --samp-prev 0.055,nan --pop-prev 0.05,nan --out chip_vs_ebmd   # 患病率数值只是示例
```

使用时的几个坑：

1. **尺度。** 二分类性状的遗传度有"观察尺度"和"liability 尺度"之分，后者需要人群患病率。[10 数据与注释](10-data-druggability-aac.md) 显示，报告用到的 CHIP 汇总统计文件里有 `CHISQ_BOLT_LMM` 一列，说明是混合线性模型（BOLT-LMM）的结果，beta 是线性尺度，不是 log OR。
2. **人群。** 欧洲人群的 LD 分数只适用于欧洲人群的 GWAS。
3. **SNP 集合。** 通常只用 HapMap3 的约 120 万个 SNP，并排除 MHC 区域。
4. **样本量。** 病例对照要用有效样本量（或分别给病例数和对照数）。

### 3.8 HDL（高清似然法）

报告里的 HDL 是 High-Definition Likelihood（Ning、Pawitan、Shen，Nat Genet 2020）。它同样用汇总统计估计遗传度和遗传相关，但不像 LDSC 那样只用 LD 分数这一个一维摘要，而是用完整的 LD 矩阵写出 z 分数的联合似然。信息利用得更充分，原文报告相比 LDSC 估计方差降低约 60%。代价是要用专门的 LD 参考面板，计算更重。

两者互为交叉验证，报告同时用两种方法是合理的。注意它们对同一份数据的结论应当方向一致，数字不必相同。

### 小结

- LDSC 把 χ² 对 LD 分数回归：斜率给遗传度，截距给混杂。
- 两个性状：斜率给遗传协方差，截距给样本重叠；所以 rg 对样本重叠稳健。
- CHIP 的有效样本量只有几万到几千，功效是这一层最大的风险。
- 单次估计噪声大，要看标准误，也要明白"没有信号"可能只是功效不足。

---

[← 上一篇：01 · CHIP 与 VAF：数据从哪来](01-chip-and-vaf.md) · [目录](README.md) · [下一篇：03 · 局部遗传相关：LAVA 与 SUPERGNOVA →](03-local-genetic-correlation.md)
