[← 上一篇：08 · 虚拟扰动：CellOracle、scTenifoldKnk、NicheNet 与 Geneformer](08-in-silico-perturbation.md) · [目录](README.md) · [下一篇：10 · 数据与注释：GWAS Catalog、PubMed、Open Targets、影像 AAC 评分 →](10-data-druggability-aac.md)

---

# 09 · 个体级队列分析：CHIP × 骨密度 × 骨折，以及共病

> 这一章带走什么：UK Biobank 里"CHIP 状态 × 骨密度 × 骨折"的分析怎么设计；年龄这个最大的混杂；按基因和克隆大小分层，样本够不够；共病分析里"观察/期望比"怎么算、怎么读，什么时候会骗你。

### 10.1 要做什么

报告第二阶段申请 UK Biobank Tier 1 数据（报告列的费用是 £3,000，学生 £500），做个体级分析：

- CHIP 状态（字段 30106、30107：变异与 VAF）× 骨密度（跟骨 BMD 约 27.9 万人；DXA 腰椎等约 6.8 万人）× 骨折结局；
- 按基因（DNMT3A、TET2、ASXL1……）和克隆大小（VAF 是否不低于 0.1）分层；
- 共病分析：携带者中"低骨密度 + 心血管事件"同时出现的比例，是否高于两个单独风险的乘积。

报告估计跟骨 BMD 人群里约有 9,000 名携带者，DXA 人群里约有 2,200 到 2,300 名。（以上数字引自报告，我没有核实。）

### 10.2 造一个玩具队列

10 万人，携带率约 3.3%（由报告里"跟骨 BMD 人群约 9,000 名携带者、27.9 万人"估出）。设计：

- CHIP 携带概率随年龄上升；
- 携带者中：DNMT3A 55%、TET2 18%、ASXL1 7%、其他 20%（常见的大致构成，仅为设定）；
- VAF 服从对数正态，使大克隆（VAF 不低于 0.1）约占 1/3（报告引用的 CHIP GWAS 里大克隆病例约占四成）；
- 骨密度（单位：标准差）随年龄下降、女性更低；CHIP 的真实效应设为：小克隆 DNMT3A −0.03、TET2 −0.08、ASXL1 −0.08、其他 0；大克隆在此基础上再减 0.12、0.18、0.20、0.05。**这些全是设定，不是真实效应大小。**

```r
set.seed(71)
N <- 100000
age <- runif(N, 40, 70); female <- rbinom(N, 1, 0.5); bmi <- rnorm(N, 27, 4); smoke <- rbinom(N, 1, 0.12)
a0 <- uniroot(function(a) mean(plogis(a + 0.12 * (age - 55))) - 0.033, c(-8, 0))$root   # 让总体携带率约 3.3%
chip <- rbinom(N, 1, plogis(a0 + 0.12 * (age - 55)))                                      # 携带概率随年龄上升
gene <- rep("无", N)
gene[chip == 1] <- sample(c("DNMT3A", "TET2", "ASXL1", "其他"), sum(chip), TRUE, prob = c(0.55, 0.18, 0.07, 0.20))
vaf   <- ifelse(chip == 1, pmin(rlnorm(N, log(0.065), 0.9), 0.45), 0)
large <- as.integer(vaf >= 0.10)                                                          # 大克隆：VAF >= 0.1
eff_small <- c(DNMT3A = -0.03, TET2 = -0.08, ASXL1 = -0.08, 其他 = 0)                     # 真实效应（骨密度的 SD 单位）
eff_large <- c(DNMT3A = -0.12, TET2 = -0.18, ASXL1 = -0.20, 其他 = -0.05)                 # 大克隆在此基础上再减
chip_eff  <- ifelse(chip == 1, eff_small[gene] + large * eff_large[gene], 0)
bmd <- 0.6 - 0.035 * (age - 55) - 0.35 * female - 0.02 * (age - 55) * female +
       0.03 * (bmi - 27) - 0.20 * smoke + chip_eff + rnorm(N)
group <- ifelse(chip == 0, "无CHIP", paste0(gene, ifelse(large == 1, "-大克隆", "-小克隆")))
group <- relevel(factor(group), ref = "无CHIP")
c(携带者 = sum(chip), 携带率 = round(mean(chip), 3), 大克隆占比 = round(mean(large[chip == 1]), 2))
```

```text
    携带者     携带率 大克隆占比 
  3381.000      0.034      0.320 
```

### 10.3 最大的混杂：年龄

CHIP 随年龄上升，骨密度也随年龄下降。所以如果不调整年龄，会把"年龄的作用"算到 CHIP 头上。

```r
fit0 <- lm(bmd ~ chip)                                                    # 不调整
fit1 <- lm(bmd ~ chip + age + female + bmi + smoke)                       # 调整年龄、性别、BMI、吸烟
ci <- function(f) round(c(估计 = coef(f)["chip"], confint(f)["chip", ]), 3)
rbind(不调整 = ci(fit0), 调整后 = ci(fit1))
cat("\n真实的平均 CHIP 效应 =", round(mean(chip_eff[chip == 1]), 3), "\n")
```

```text
       估计.chip  2.5 % 97.5 %
不调整    -0.442 -0.480 -0.404
调整后    -0.101 -0.136 -0.066

真实的平均 CHIP 效应 = -0.075 
```

| 模型 | CHIP 的估计效应（SD） | 95% 区间 |
|---|---|---|
| 不调整 | −0.442 | −0.480 到 −0.404 |
| 调整年龄、性别、BMI、吸烟 | −0.101 | −0.136 到 −0.066 |
| 真值（携带者的平均） | −0.075 | |

不调整时，估计出的效应是真值的近 6 倍。调整后，区间覆盖真值。

调整变量的原则：只调整**不在因果路径上的混杂**（年龄、性别、BMI、吸烟、评估中心、遗传主成分、绝经和激素治疗、糖皮质激素使用等），不要调整"CHIP 的后果"，如已经确诊的骨质疏松或已在用的抗骨质疏松药。

### 10.4 按基因和克隆大小分层

```r
fit2 <- lm(bmd ~ group + age + female + bmi + smoke)
tab  <- summary(fit2)$coefficients; keep <- grep("^group", rownames(tab))
res  <- data.frame(估计 = tab[keep, 1], SE = tab[keep, 2], row.names = sub("^group", "", rownames(tab)[keep]))
res$真实效应 <- sapply(rownames(res), function(g) { gg <- strsplit(g, "-")[[1]]; eff_small[gg[1]] + (gg[2] == "大克隆") * eff_large[gg[1]] })
res$携带者数 <- table(group)[rownames(res)]
round(res[order(rownames(res)), ], 3)
```

```text
                估计    SE 真实效应 携带者数
ASXL1-大克隆  -0.103 0.122    -0.28       68
ASXL1-小克隆  -0.131 0.079    -0.08      163
DNMT3A-大克隆 -0.209 0.042    -0.15      580
DNMT3A-小克隆 -0.065 0.028    -0.03     1311
TET2-大克隆   -0.269 0.071    -0.26      203
TET2-小克隆   -0.018 0.051    -0.08      396
其他-大克隆   -0.079 0.068    -0.05      222
其他-小克隆   -0.065 0.048     0.00      438
```

读法：

- 携带者越多，估计越准：DNMT3A 小克隆（1,311 人）SE 0.028，ASXL1 大克隆（68 人）SE 0.122；
- 大克隆的估计普遍更负，与设定方向一致；
- 小亚组的估计会偏离真值：ASXL1 小克隆估计 −0.131（真值 −0.08），TET2 小克隆 −0.018（真值 −0.08），都落在一个 SE 的范围内。

**两个常见的误判：**

1. "一个亚组显著、另一个不显著，所以两者有差别"是错的。两个亚组的差别要用**交互项或对比检验**直接检验（Gelman 与 Stern，Am Stat 2006）。
2. 8 个亚组就是 8 次检验，要做多重检验校正（FDR 或 Bonferroni），并事先写好哪些是主要亚组。

### 10.5 样本够不够：最小可检出差异

假设携带率 3.3%，携带者构成如下（假设值），骨密度标准差为 1，80% 功效、双侧 5%：

```r
# 在 27.9 万（跟骨 BMD）和 6.8 万（DXA）两种规模下，各亚组能检出多大的差异（骨密度 SD = 1，80% 功效，双侧 5%）
share <- c(DNMT3A = 0.55, TET2 = 0.18, ASXL1 = 0.07, `大克隆（所有基因）` = 0.30)     # 假设的携带者构成
prev  <- 0.033
mdd <- function(Ntot, s) { nc <- Ntot * prev * s; 2.8 * sqrt(1 / nc + 1 / (Ntot * (1 - prev))) }
tab <- data.frame(跟骨BMD_携带者数 = round(279000 * prev * share), 跟骨BMD_最小可检出差异 = round(mdd(279000, share), 2),
                  DXA_携带者数 = round(68000 * prev * share),      DXA_最小可检出差异 = round(mdd(68000, share), 2),
                  row.names = names(share))
tab
```

```text
                   跟骨BMD_携带者数 跟骨BMD_最小可检出差异 DXA_携带者数 DXA_最小可检出差异
DNMT3A                         5064                   0.04         1234               0.08
TET2                           1657                   0.07          404               0.14
ASXL1                           644                   0.11          157               0.22
大克隆（所有基因）             2762                   0.05          673               0.11
```

读法：

- 跟骨 BMD（27.9 万人）：DNMT3A 能检出 0.04 个 SD 的差异，TET2 0.07，ASXL1 0.11；
- DXA（6.8 万人）：DNMT3A 0.08，TET2 0.14，ASXL1 0.22。

报告说"9,000 名携带者足以做按基因分层的分析"，在跟骨 BMD 上基本成立；但 DXA 人群里 ASXL1 只有约 157 人，只能检出 0.22 个 SD 的差异，要做好"可能没有结果"的心理准备。DXA 是金标准但样本小，跟骨超声样本大但是估计值，两者要互相对照。

### 10.6 骨折：Cox 比例风险模型

骨折是"时间到事件"的结局，用 Cox 模型。设定：骨折风险随年龄上升、女性更高、骨密度越低越高；大克隆另有一个直接效应。随访 10 年。

```r
library(survival)
lin  <- 0.07 * (age - 55) + 0.35 * female - 0.30 * bmd + 0.15 * large                 # 骨折风险：低骨密度升高风险，大克隆另有直接作用
t_ev <- rexp(N, 0.004 * exp(lin))                                                    # 事件时间（年）
time <- pmin(t_ev, 10); status <- as.integer(t_ev <= 10)                             # 随访 10 年
chip_lv <- factor(ifelse(chip == 0, "无CHIP", ifelse(large == 1, "大克隆", "小克隆")), levels = c("无CHIP", "小克隆", "大克隆"))
fit <- coxph(Surv(time, status) ~ chip_lv + age + female + bmi + smoke)
c(骨折例数 = sum(status))
round(summary(fit)$conf.int[1:2, c(1, 3, 4)], 2)
```

```text
骨折例数 
    5696 
              exp(coef) lower .95 upper .95
chip_lv小克隆      1.04      0.91      1.19
chip_lv大克隆      1.24      1.04      1.49
```

大克隆的风险比 1.24（1.04 到 1.49），小克隆 1.04（0.91 到 1.19）。

使用要点：

1. **随访起点：** 从取血或外显子测序之日起算，不要把之前发生的骨折算进来；
2. **比例风险假设：** 要检验（如 `cox.zph`）；
3. **竞争风险：** 死亡会阻止骨折发生，高龄人群里要考虑 Fine-Gray 模型；
4. **是否调整骨密度：** 调整前的 HR 是 CHIP 的"总效应"；把骨密度也放进模型，HR 会减弱，减弱的部分对应"通过骨密度起作用"的路径，这是中介分析的问题，要事先想清楚。

### 10.7 共病分析：观察/期望比

报告要问："CHIP 携带者中'低骨密度 + 心血管事件'的联合发生率，是否高于两个单独风险的乘积？"

**为什么看乘积。** 如果两件事互相独立，同时发生的概率就等于各自概率的乘积。同时出现的次数比这个"期望"多，说明两件事在这批人里**结伴出现**。

**怎么算。** 观察/期望比 = 观察到的同时发生人数 ÷ 期望人数。期望 = 每个分层里 人数 × P(低骨密度) × P(心血管事件)，再加总。

**为什么要分层。** 年龄会让两件事天然结伴出现（越老越容易两者都有）。所以要在**年龄 × 性别**的分层里算，把年龄的作用拿掉。

**怎么读。** 只看携带者的比值是否大于 1 不够，还要和**非携带者**比。因为吸烟、肾功能、衰弱等共同原因，在没有 CHIP 的人群里也会让两件事结伴出现。

玩具模拟：两个世界，边际风险完全相同，区别只在**携带者内部的依赖结构**：

- 共享同一套程序：携带者体内有一个"程序强度"，同时推高低骨密度和心血管事件的风险；
- 各自独立：两件事各自有独立的"强度"，边际风险一样。

```r
# 共病：低骨密度和心血管事件是否“扎堆”出现在携带者身上？
sim_comorb <- function(shared) {
  u1 <- rnorm(N); u2 <- if (shared) u1 else rnorm(N)     # 携带者体内的“程序强度”：共享同一个 / 各自独立
  low <- rbinom(N, 1, plogis(-1.1 + 0.05 * (age - 55) + 0.5 * female + 0.3 * smoke + chip * (0.35 + 0.8 * u1)))
  cvd <- rbinom(N, 1, plogis(-2.8 + 0.07 * (age - 55) + 0.3 * (1 - female) + 0.5 * smoke + chip * (0.25 + 0.8 * u2)))
  data.frame(age, female, chip, low, cvd)
}
oe <- function(d) {                                       # 在“年龄 x 性别”分层里：观察到的同时出现 vs 各自风险相乘的期望
  st <- interaction(cut(d$age, seq(40, 70, 5), include.lowest = TRUE), d$female)
  O <- sum(d$low & d$cvd)
  E <- sum(sapply(split(d, st), function(s) nrow(s) * mean(s$low) * mean(s$cvd)))
  c(观察 = O, 期望 = E, 比值 = O / E)
}
boot_ci <- function(d, B = 300) quantile(replicate(B, oe(d[sample(nrow(d), replace = TRUE), ])[3]), c(0.025, 0.975))
summarize <- function(d) {
  car <- oe(d[d$chip == 1, ]); non <- oe(d[d$chip == 0, ]); ci <- boot_ci(d[d$chip == 1, ])
  c(携带者_观察 = car[[1]], 携带者_期望 = round(car[[2]], 1), 携带者_比值 = round(car[[3]], 2),
    比值_下限 = round(ci[[1]], 2), 比值_上限 = round(ci[[2]], 2), 非携带者_比值 = round(non[[3]], 2))
}
set.seed(72)
rbind(`共享同一套程序` = summarize(sim_comorb(TRUE)), `各自独立` = summarize(sim_comorb(FALSE)))
```

```text
               携带者_观察 携带者_期望 携带者_比值 比值_下限 比值_上限 非携带者_比值
共享同一套程序         336       274.4        1.22      1.15      1.29          1.03
各自独立               298       297.5        1.00      0.93      1.07          0.99
```

| 世界 | 携带者 观察 / 期望 | 携带者 O/E（95% 区间） | 非携带者 O/E |
|---|---|---|---|
| 共享同一套程序 | 336 / 274.4 | 1.22（1.15 到 1.29） | 1.03 |
| 各自独立 | 298 / 297.5 | 1.00（0.93 到 1.07） | 0.99 |

读法：共享世界里，携带者的 O/E 是 1.22，明显大于 1，而非携带者只有 1.03；独立世界里携带者的 O/E 是 1.00。两个世界里携带者的边际风险一样，但 O/E 能把它们区分开。

**但 O/E > 1 不能证明"同一套破骨样程序"。** 任何让两件事结伴的共同原因都会造成 O/E > 1（炎症、肾功能、衰弱、药物……）。要加强这条证据，可以：

1. 与非携带者比较（做法如上）；
2. 看剂量反应：大克隆的 O/E 是否高于小克隆；
3. 尽可能调整更多的共同原因；
4. 更严格地，用多状态模型或多项逻辑回归，联合建模"都没有、只有骨、只有血管、两者都有"四种结局。

### 10.8 UKB-RAP 的限制

报告提醒：所有分析都在 UKB-RAP 云平台上进行，计算与存储按量计费；数据不能下载到本地。所以代码要写成可复现、可批量提交的形式，先在小样本上调试，再放到完整数据上跑。

### 小结

- 不调整年龄，CHIP 对骨密度的效应会被放大近 6 倍。
- 分层要比较"交互"，不是比较"谁显著"；亚组多了要做多重检验。
- 跟骨 BMD 的样本基本够，DXA 里 ASXL1 只能检出大效应。
- O/E 比必须在年龄 × 性别分层里算，并与非携带者比较；O/E > 1 也不能单独证明共享程序。

---

[← 上一篇：08 · 虚拟扰动：CellOracle、scTenifoldKnk、NicheNet 与 Geneformer](08-in-silico-perturbation.md) · [目录](README.md) · [下一篇：10 · 数据与注释：GWAS Catalog、PubMed、Open Targets、影像 AAC 评分 →](10-data-druggability-aac.md)
