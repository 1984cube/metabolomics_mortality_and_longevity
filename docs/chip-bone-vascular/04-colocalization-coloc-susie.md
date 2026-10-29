[← 上一篇：03 · 局部遗传相关：LAVA 与 SUPERGNOVA](03-local-genetic-correlation.md) · [目录](README.md) · [下一篇：05 · 单细胞分析基础，以及“别把细胞当重复” →](05-single-cell-basics-pseudobulk.md)

---

# 04 · 共定位：coloc 与 SuSiE-coloc

> 这一章带走什么：为什么"两个性状在同一个区域都显著"不等于"同一个变异"；coloc 怎么算五个后验概率；单一因果变异的假设什么时候会骗你；SuSiE-coloc 怎么补救。

### 5.1 要回答的问题

一个基因组区域里，一个因果变异会让它周围一串 SNP 都显著（因为 LD）。设想 CHIP 和骨密度在同一个区域都有显著信号：

- 可能是**同一个变异**同时影响两个性状（这是我们想要的）；
- 也可能是**两个不同的变异**碰巧靠得很近，各影响一个性状。

只看"都显著"分不清这两种情形。共定位（colocalization）就是来判断：两个性状的信号来自同一个因果变异的概率有多大。

### 5.2 五个假设

对一个区域，coloc 比较五个互斥的假设：

| 假设 | 含义 | 默认先验 |
|---|---|---|
| H0 | 两个性状都没有关联 | 其余的补足 |
| H1 | 只有性状 1 有关联 | 每个 SNP 对性状 1 有因果的概率 p1 = 1e-4 |
| H2 | 只有性状 2 有关联 | p2 = 1e-4 |
| H3 | 两个性状都有关联，但因果变异不同 | p1 × p2 |
| H4 | 两个性状都有关联，共享同一个因果变异 | p12 = 1e-5 |

注意 H4 与 H3 的先验比：对每一对 SNP，H4 的先验是 1e-5，H3 的先验是 1e-4 × 1e-4 = 1e-8，**H4 比 H3 有 1000 倍的先验优势**。这个优势会在后面的"多信号"例子里起作用。

输出是五个后验概率 PP.H0 到 PP.H4，加起来等于 1。

### 5.3 近似贝叶斯因子

coloc 的核心是给每个 SNP 算一个"贝叶斯因子"：这个 SNP 有因果效应的证据有多强。用 Wakefield 的近似，只需要 beta 和 SE：

```text
V = SE²                              （效应估计的方差）
W = (0.15 × sdY)²                    （先验：真实效应的方差，连续性状的默认值）
r = W ÷ (W + V)
z = beta ÷ SE
log ABF = 0.5 × ( log(1 − r) + r × z² )
```

手算一个例子：beta = 0.25，SE = 0.025，sdY = 1。

```text
V = 0.025² = 0.000625；  W = 0.15² = 0.0225
r = 0.0225 ÷ (0.0225 + 0.000625) = 0.973
z = 0.25 ÷ 0.025 = 10
log ABF = 0.5 × ( log(0.027) + 0.973 × 100 ) = 0.5 × ( −3.61 + 97.30 ) = 46.84
```

z 越大，log ABF 越大，证据越强。

### 5.4 汇总成五个后验

对每个性状，把区域内所有 SNP 的 ABF 加起来得到 H1、H2 的证据；把两个性状同一 SNP 的 ABF 相乘再加起来得到 H4 的证据；"所有 SNP 对"的乘积之和减去"同一 SNP"的部分得到 H3 的证据。再乘上先验，归一化。这些都在对数尺度上算（避免数值溢出）。下面是完整的实现，不到 20 行：

```r
# coloc 的核心：每个 SNP 的近似贝叶斯因子（Wakefield），再汇总成五个假设的后验概率
log_abf <- function(beta, varbeta, sdY = 1) {
  W <- (0.15 * sdY)^2                                 # 先验：效应大小的标准差 = 0.15 个 sdY
  r <- W / (W + varbeta); z <- beta / sqrt(varbeta)
  0.5 * (log(1 - r) + r * z^2)
}
logsum  <- function(x) { m <- max(x); m + log(sum(exp(x - m))) }
logdiff <- function(x, y) { m <- max(x, y); m + log(exp(x - m) - exp(y - m)) }
my_coloc <- function(d1, d2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5) {
  l1 <- log_abf(d1$beta, d1$varbeta, d1$sdY); l2 <- log_abf(d2$beta, d2$varbeta, d2$sdY)
  lH <- c(0,                                                         # H0：都没有关联
          log(p1) + logsum(l1),                                      # H1：只有性状 1
          log(p2) + logsum(l2),                                      # H2：只有性状 2
          log(p1) + log(p2) + logdiff(logsum(l1) + logsum(l2), logsum(l1 + l2)),   # H3：各有因果变异，但不是同一个
          log(p12) + logsum(l1 + l2))                                # H4：共享同一个因果变异
  setNames(exp(lH - logsum(lH)), paste0("H", 0:4))
}
```

### 5.5 实战：三种情形

造一个区域：150 个 SNP，相邻 SNP 相关 0.95，向两侧衰减（LD 很强）。每个性状来自 4000 人的独立样本，所以是"两样本"设计。

```r
set.seed(7)
m <- 150; n <- 4000                                   # 一个区域 150 个 SNP；每个研究 4000 人
chol_R <- chol(0.95^abs(outer(1:m, 1:m, "-")))        # 相邻 SNP 相关 0.95，向两侧衰减（LD 很强）
f   <- runif(m, 0.15, 0.50)
snp <- sprintf("s%03d", 1:m)
sim_geno <- function(n) {                             # 造这个区域的基因型 0/1/2
  Z <- matrix(rnorm(n * m), n, m) %*% chol_R
  (Z > rep(qnorm((1 - f)^2), each = n)) + (Z > rep(qnorm(1 - f^2), each = n))
}
make_study <- function(causal, beta = 0.25) {         # 造一个研究：因果 SNP 的位置 + 逐 SNP 回归
  G <- sim_geno(n)
  y <- rnorm(n) + if (length(causal)) as.vector(G[, causal, drop = FALSE] %*% rep_len(beta, length(causal))) else 0
  y <- as.vector(scale(y))                            # 标准化，sdY = 1
  gc <- sweep(G, 2, colMeans(G)); sxx <- colSums(gc^2); sxy <- as.vector(crossprod(gc, y))
  b <- sxy / sxx; se <- sqrt((sum(y^2) - b * sxy) / (n - 2) / sxx)
  list(beta = b, varbeta = se^2, snp = snp, type = "quant", sdY = 1, N = n)
}
Gref <- sim_geno(2000); LD <- cor(Gref); dimnames(LD) <- list(snp, snp)   # LD 参考样本
```

三种情形：

- **H4（共享）：** 两个性状的因果 SNP 都在位置 60；
- **H3（不同）：** 性状 1 的因果 SNP 在 60，性状 2 的在 80（两者的 LD 中等）；
- **H1（只有性状 1）：** 性状 2 没有因果 SNP。

```r
d_shared <- list(make_study(60), make_study(60))      # H4：两个性状共享因果 SNP（位置 60）
d_diff   <- list(make_study(60), make_study(80))      # H3：因果 SNP 不同（60 和 80，LD 中等）
d_only1  <- list(make_study(60), make_study(NULL))    # H1：只有性状 1 有因果 SNP
cases <- list("共享(H4)" = d_shared, "不同(H3)" = d_diff, "只有性状1(H1)" = d_only1)
round(t(sapply(cases, function(d) my_coloc(d[[1]], d[[2]]))), 3)
```

```text
              H0    H1 H2    H3    H4
共享(H4)       0 0.000  0 0.000 1.000
不同(H3)       0 0.000  0 1.000 0.000
只有性状1(H1)  0 0.962  0 0.005 0.034
```

三种情形都判对了：共享时 PP.H4 = 1.000；不同时 PP.H3 = 1.000；只有性状 1 时 PP.H1 = 0.962。

和软件包对账：同样的输入，我手写的函数与 `coloc` 包的 `coloc.abf` 结果应当完全一致。

```r
# 和软件包对账：同样的输入，结果应当完全一致
library(coloc)
pkg_pp <- function(d) {
  invisible(capture.output(res <- suppressMessages(coloc.abf(dataset1 = d[[1]], dataset2 = d[[2]]))))
  unname(res$summary[paste0("PP.H", 0:4, ".abf")])
}
res_pkg <- t(sapply(cases, pkg_pp))
dimnames(res_pkg) <- list(names(cases), paste0("H", 0:4))
round(res_pkg, 3)
all.equal(unname(res_pkg), unname(t(sapply(cases, function(d) my_coloc(d[[1]], d[[2]])))))
```

```text
              H0    H1 H2    H3    H4
共享(H4)       0 0.000  0 0.000 1.000
不同(H3)       0 0.000  0 1.000 0.000
只有性状1(H1)  0 0.962  0 0.005 0.034
[1] TRUE
```

`all.equal` 返回 TRUE，说明手写的实现与软件包一致。

看一眼数据，了解"信号有多强、峰在哪"：

```r
# 看一眼数据：每个性状最强的 SNP 在哪、z 有多大
top <- function(d) c(位置 = which.max(abs(d$beta)), z = round(max(abs(d$beta / sqrt(d$varbeta))), 1))
rbind(共享_性状1 = top(d_shared[[1]]), 共享_性状2 = top(d_shared[[2]]),
      不同_性状1 = top(d_diff[[1]]),   不同_性状2 = top(d_diff[[2]]))
```

```text
           位置    z
共享_性状1   60  9.5
共享_性状2   60  8.9
不同_性状1   60 11.1
不同_性状2   80 10.9
```

"不同"情形里，两个性状的最强信号分别在位置 60 和 80，z 都超过 10，一眼就能看出峰不在同一个地方。

### 5.6 单一因果变异的假设会骗你

coloc 默认**每个性状在这个区域里最多一个因果变异**。如果一个性状有两个独立信号，会怎样？

造这样一个情形：性状 1 有两个因果 SNP：位置 60（强）和位置 100（弱）；性状 2 只有位置 100。也就是说，**真正共享的是性状 1 的"次要信号"**。

```r
# 单一因果变异假设被打破：性状 1 有两个因果 SNP（60 强、100 弱），性状 2 只有 100
# 也就是：真正共享的是性状 1 的“次要信号”
d_multi <- list(make_study(c(60, 100), beta = c(0.35, 0.12)), make_study(100))
rbind(性状1 = top(d_multi[[1]]), 性状2 = top(d_multi[[2]]))
round(my_coloc(d_multi[[1]], d_multi[[2]]), 3)
```

```text
      位置    z
性状1   60 15.2
性状2  100 11.1
H0 H1 H2 H3 H4 
 0  0  0  1  0 
```

性状 1 的最强信号在 60（z = 15.2），性状 2 的在 100（z = 11.1）。coloc 看到两个性状的"主峰"不在同一处，判成 H3（因果变异不同，PP.H3 = 1）。**这是错的**：它们在位置 100 上确实共享，只是被性状 1 的强信号盖住了。

为什么 H4 没有赢？H4 的证据来自"同一个 SNP 上两个性状都强"，只能用位置 100。但性状 1 在位置 100 的证据（z = 6.5，log ABF = 18.9）远不如位置 60（z = 15.2，log ABF = 111.1），差了 92；而 H4 相对 H3 的 1000 倍先验优势，在对数尺度上只有 ln(1000) = 6.9。差距太大，先验补不回来。

### 5.7 SuSiE-coloc

SuSiE（Sum of Single Effects）是一种精细定位方法：把一个区域里的信号拆成若干个"单一效应"之和，每个单一效应对应一个**可信集**（一小撮 SNP，包含因果变异的概率不低于 95%）。SuSiE-coloc（Wallace 2021，PLoS Genet）先用 SuSiE 分别给两个性状找出可信集，再对**每一对可信集**做共定位。

它需要一个 LD 矩阵（来自与 GWAS 相同人群的参考样本）。下面用前面造的 LD 参考样本。

```r
library(susieR)
d1 <- d_multi[[1]]; d2 <- d_multi[[2]]; d1$LD <- LD; d2$LD <- LD
S1 <- invisible(suppressMessages(runsusie(d1))); S2 <- invisible(suppressMessages(runsusie(d2)))
c(性状1的可信集个数 = length(S1$sets$cs), 性状2的可信集个数 = length(S2$sets$cs))
res <- as.data.frame(coloc.susie(S1, S2)$summary)
res[, c("hit1", "hit2", "PP.H3.abf", "PP.H4.abf")]
```

```text
性状1的可信集个数 性状2的可信集个数 
                2                 1 
  hit1 hit2    PP.H3.abf    PP.H4.abf
1 s060 s100 1.0000000000 3.810621e-22
2 s100 s100 0.0004530253 9.994414e-01
```

读法：

- 性状 1 找到 2 个可信集（对应位置 60 和 100），性状 2 找到 1 个；
- 配对 (s060, s100)：PP.H3 ≈ 1，不是同一个变异；
- 配对 (s100, s100)：**PP.H4 = 0.9994**，共享同一个变异。

SuSiE-coloc 正确找到了被主信号盖住的共享变异。

### 5.8 实战要点

1. **等位基因要对齐。** 两个性状的汇总统计必须是同一批 SNP、同一个效应等位基因（见 [MR 讲义](../mendelian-randomization.md) 第 5 章）。
2. **区域选择。** 通常取每个显著位点附近（如两侧 500 kb）的所有 SNP。两个性状的区域坐标版本要一致（GRCh37 或 GRCh38）。
3. **阈值。** 常用经验线是 PP.H4 不低于 0.8；也可以同时看 PP.H3 + PP.H4（有没有信号）和 PP.H4 ÷ PP.H3。PP.H0、PP.H1、PP.H2 高时，往往只是功效不足，不是"没有共定位"。
4. **先验敏感性。** 把 p12 在 1e-6 到 1e-4 之间变化，看结论是否稳定（`coloc` 包有 `sensitivity()` 函数）。
5. **二分类性状。** CHIP 这类病例对照性状，如果汇总统计来自 BOLT-LMM 线性混合模型，beta 是线性尺度，可以按连续性状处理，`sdY` 取 √(p(1−p))，p 是病例比例。具体做法以 `coloc` 文档为准。
6. **三个以上性状。** 报告要找"三表型共享位点"。可以两两做共定位再取交集；多性状共定位有专门的方法（如 HyPrColoc，Foley 等，Nat Commun 2021）。
7. **LD 必须匹配。** SuSiE-coloc 的 LD 矩阵与 GWAS 人群不匹配时，会找出假的可信集。

### 小结

- 共定位回答"同一个变异"的问题，输出 PP.H0 到 PP.H4。
- 核心是近似贝叶斯因子 + 五个假设的先验；手写不到 20 行，与 `coloc` 包完全一致。
- 单一因果变异的假设在多信号区域会骗你；SuSiE-coloc 按可信集配对，能找回被盖住的共享信号。

---

[← 上一篇：03 · 局部遗传相关：LAVA 与 SUPERGNOVA](03-local-genetic-correlation.md) · [目录](README.md) · [下一篇：05 · 单细胞分析基础，以及“别把细胞当重复” →](05-single-cell-basics-pseudobulk.md)
