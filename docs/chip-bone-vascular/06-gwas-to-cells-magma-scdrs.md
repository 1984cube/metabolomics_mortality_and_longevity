[← 上一篇：05 · 单细胞分析基础，以及“别把细胞当重复”](05-single-cell-basics-pseudobulk.md) · [目录](README.md) · [下一篇：07 · 基因程序评分与跨组织比对 →](07-gene-program-scoring.md)

---

# 06 · 把 GWAS 信号映射到细胞：MAGMA、scDRS 与 scPagwas

> 这一章带走什么：怎么把"基因组上的关联信号"和"单细胞图谱里的表达"接起来；scDRS 为什么必须用"匹配对照"；用一个玩具例子看匹配和不匹配的差别。

### 7.1 要回答的问题

GWAS 告诉你：这个性状与基因组上很多位点有关。单细胞图谱告诉你：每个细胞表达哪些基因。把两者接起来，可以问：**这个性状的遗传风险，偏向哪一类细胞？**

对应报告第一阶段的第 4 步：把 AAC、CAC、eBMD 的 GWAS 信号映射到钙化斑块和骨髓单细胞图谱，看风险是否富集在髓系细胞上。

### 7.2 第一步：从 SNP 到基因（MAGMA）

GWAS 的单位是 SNP，单细胞的单位是基因，要先对齐。MAGMA（de Leeuw 等，PLoS Comput Biol 2015）做基因水平分析：把一个基因及其附近 SNP 的 p 值合并成一个**基因水平的 z 分数**，同时考虑 LD 和基因大小。z 越大，说明这个基因越可能与性状有关。

### 7.3 第二步：给每个细胞打"疾病分数"（scDRS）

scDRS（Zhang 等，Nat Genet 2022）的思路：

1. 取 GWAS 排名靠前的一批基因（论文默认约 1,000 个），用 MAGMA 的 z 分数作权重；
2. 每个细胞的"原始疾病分数"= 这些基因表达量的加权平均；
3. 但这个分数不能直接比。要生成很多组（论文默认约 1,000 组）**表达量与方差匹配**的随机对照基因集，同样算分数，作为"没有疾病信号时的期望"；
4. 用对照分数把疾病分数标准化，得到每个细胞的标准化分数和 p 值；
5. 在细胞类型层面汇总，如"某类细胞里高分细胞的比例"，并看类型内部的异质性。

以上是论文设计的概述，默认参数和细节以原文及软件文档为准。

### 7.4 为什么对照必须"匹配"

GWAS 筛出来的基因，往往本来表达量就偏高（高表达基因更容易被检出关联，也更常出现在功能基因里）。如果随机对照基因集是从**所有**基因里抽的，它们大部分是低表达基因，疾病基因集的分数永远比对照高，**所有细胞类型都会"显著"**。所以对照基因集必须在表达量上与疾病基因集匹配。

### 7.5 实战：玩具 scDRS

造一个玩具图谱：4 种细胞类型（T 细胞、B 细胞、髓系、基质）各 1,000 个细胞，2,000 个基因。设计：

- 疾病基因集：100 个基因，都来自表达量最高的前 400 个基因（模拟真实 GWAS 基因偏高表达）；
- 真信号：其中 60 个基因在髓系细胞里表达更高；
- 对照性状：另外 100 个同样高表达的基因，但没有任何细胞类型特异性。

```r
set.seed(41)
ngene <- 2000; per <- 1000; types <- c("T细胞", "B细胞", "髓系", "基质")
type  <- rep(types, each = per); ncell <- length(type)
mu_g  <- rnorm(ngene, -0.5, 1.4)                         # 每个基因的基础表达量（对数尺度），差别很大
fc    <- matrix(0, ngene, 4)                             # 每种细胞类型里 300 个基因有特异的上调或下调
for (k in 1:4) fc[sample(ngene, 300), k] <- rnorm(300, 0, 0.8)

hi  <- order(mu_g, decreasing = TRUE)[1:400]             # 表达量最高的 400 个基因
dis <- sample(hi, 100)                                   # 疾病基因集：100 个，都来自高表达基因（真实 GWAS 基因常如此）
sig <- sample(dis, 60); fc[sig, 3] <- fc[sig, 3] + 0.8   # 真信号：其中 60 个在髓系细胞里更高
dis_null <- sample(setdiff(hi, dis), 100)                # 对照性状：同样是 100 个高表达基因，但没有任何细胞类型特异性
w <- runif(100, 3, 8)                                    # 基因权重（真实分析里来自 MAGMA 基因水平 z 分数）

depth  <- rlnorm(ncell, 0, 0.3)
lambda <- depth * exp(matrix(mu_g, ncell, ngene, byrow = TRUE) + t(fc)[match(type, types), ])
counts <- matrix(rpois(length(lambda), lambda), ncell)
X <- log1p(counts / rowSums(counts) * 1e4)               # 标准化后的表达量
```

再写打分函数：把基因按平均表达量分成 20 档；"匹配对照"是为每个疾病基因配一个同档位的随机基因；"朴素对照"是从所有基因里随便抽。

```r
raw_score <- function(genes, wt) as.vector(X[, genes, drop = FALSE] %*% (wt / sum(wt)))   # 加权平均表达量

gene_mean <- colMeans(X)
bin <- cut(gene_mean, quantile(gene_mean, 0:20 / 20), include.lowest = TRUE, labels = FALSE)   # 按平均表达量分 20 档

make_ctrl <- function(genes, wt, matched, B = 300) {     # 造 B 组对照基因集
  pool <- setdiff(seq_len(ngene), genes)
  t(sapply(seq_len(B), function(i) {
    g <- if (matched) {                                  # 匹配对照：每个疾病基因配一个同一表达档位的随机基因
      sapply(genes, function(x) sample(pool[bin[pool] == bin[x]], 1))
    } else sample(pool, length(genes))                   # 朴素对照：从所有基因里随便抽
    raw_score(g, sample(wt))
  }))
}

score_trait <- function(genes, wt, matched) {
  s    <- raw_score(genes, wt)
  ctrl <- make_ctrl(genes, wt, matched)
  z    <- (s - colMeans(ctrl)) / apply(ctrl, 2, sd)      # 每个细胞：疾病分数相对对照分数的 z
  p    <- (1 + rowSums(t(ctrl) >= s)) / (nrow(ctrl) + 1) # 蒙特卡洛 p 值：对照分数不低于疾病分数的比例
  data.frame(mean_z = tapply(z, type, mean)[types], frac_sig = tapply(p < 0.05, type, mean)[types])
}
```

对两个性状（真信号、无关性状）各用两种对照打分：

```r
set.seed(42)
r_dis_m  <- score_trait(dis, w, TRUE);        r_dis_n  <- score_trait(dis, w, FALSE)        # 真信号：匹配对照 / 朴素对照
r_null_m <- score_trait(dis_null, w, TRUE);   r_null_n <- score_trait(dis_null, w, FALSE)   # 无关性状：匹配对照 / 朴素对照
z_tab <- cbind(真信号_匹配对照 = r_dis_m$mean_z, 真信号_朴素对照 = r_dis_n$mean_z,
               无关性状_匹配对照 = r_null_m$mean_z, 无关性状_朴素对照 = r_null_n$mean_z)
p_tab <- cbind(真信号 = r_dis_m$frac_sig, 无关性状 = r_null_m$frac_sig)
rownames(z_tab) <- rownames(p_tab) <- types
cat("每种细胞类型的平均 z 分数\n"); round(z_tab, 1)
cat("\n匹配对照下，蒙特卡洛 p<0.05 的细胞比例\n"); round(p_tab, 2)
```

```text
每种细胞类型的平均 z 分数
      真信号_匹配对照 真信号_朴素对照 无关性状_匹配对照 无关性状_朴素对照
T细胞             0.5            15.7               0.4              15.3
B细胞            -0.3            14.7               1.2              16.2
髓系              5.8            19.3              -0.9              14.8
基质              0.5            15.0               0.8              15.3

匹配对照下，蒙特卡洛 p<0.05 的细胞比例
      真信号 无关性状
T细胞   0.10     0.09
B细胞   0.01     0.33
髓系    1.00     0.00
基质    0.09     0.17
```

读法：

- **匹配对照：** 真信号性状里，只有髓系细胞的平均 z 分数很高（5.8），而且 100% 的髓系细胞 p<0.05，其他三类都在 0.5 左右。无关性状里，髓系的 z 是 −0.9、没有显著细胞。方法有分辨力。
- **朴素对照：** 两个性状里，**所有**细胞类型的 z 分数都在 15 左右，所有细胞 p<0.05。这等于什么都没说。
- 提醒：无关性状里 B 细胞的 z = 1.2、33% 的细胞 p<0.05。这是随机抽的 100 个基因恰好偏向 B 细胞，属于偶然的小偏差。真实分析要在**细胞类型层面做置换检验**，并看效应大小，而不是只看显著细胞的比例。

**本玩具与真实 scDRS 的差别：** 没有对基因做标准化和协变量回归；对照基因集只有 300 组；只匹配了平均表达量，没有匹配方差；p 值逐细胞算，没有汇总各细胞的对照分数；也没有细胞类型层面的检验和异质性分析。

### 7.6 scPagwas

scPagwas（Ma 等，Cell Genomics 2023）是另一种把 GWAS 与单细胞接起来的方法，思路是在通路层面把 GWAS 信号与单细胞的通路活性联系起来。细节以原文为准，本讲义没有运行。报告写"scDRS 或 scPagwas"，建议两个都做，互相验证。

### 7.7 使用要点

1. **GWAS 的信号要够。** 信号太弱的 GWAS（如 TET2-CH，只有 2,041 例）很难给出稳定的"疾病基因集"。优先用大样本的表型（eBMD、CAC、AAC）。
2. **图谱里要有目标细胞类型。** 报告指出，破骨细胞专门的单细胞数据较少，骨这一侧是弱点。如果图谱里没有破骨样细胞，风险就不可能"富集"在它身上。
3. **富集不等于因果。** 它只说"风险基因在这类细胞里表达偏高"，不说这类细胞是原因。
4. **对应立项判据 2：** "至少一个钙化相关 GWAS 的风险富集在髓系细胞上"。要先在方案里写好什么算"富集"（细胞类型层面的 p 值阈值、校正方式）。

### 小结

- 先用 MAGMA 把 SNP 信号汇总到基因，再用 scDRS 给每个细胞打分。
- 对照基因集必须在表达量上匹配；朴素对照会让所有细胞类型都"显著"。
- 富集只说明风险偏向哪类细胞，不说明因果。

---

[← 上一篇：05 · 单细胞分析基础，以及“别把细胞当重复”](05-single-cell-basics-pseudobulk.md) · [目录](README.md) · [下一篇：07 · 基因程序评分与跨组织比对 →](07-gene-program-scoring.md)
