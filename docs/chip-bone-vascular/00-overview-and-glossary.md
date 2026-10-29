[目录](README.md) · [下一篇：01 · CHIP 与 VAF：数据从哪来 →](01-chip-and-vaf.md)

---

# CHIP × 骨–血管钙化课题：生信技术讲义（MR 除外，通俗版）

> 这是《CHIP × 骨–血管钙化悖论：课题可行性核实报告》的配套讲义。报告里除孟德尔随机化（MR）以外的生信技术，都按 MR 讲义（`docs/mendelian-randomization.md`）的写法讲：先讲直觉，再代入具体数字，最后才给公式；每一章配一个可以动手跑的 R 例子。
>
> 讲义里的 R 代码和输出由工具自动嵌入：代码在 R 4.3.3 里实际运行过，输出就是那次运行的结果，没有手抄。[04 共定位](04-colocalization-coloc-susie.md) 用到的 `coloc` 5.2.3 和 `susieR` 0.14.2 也是实际装上跑的。[10 数据与注释](10-data-druggability-aac.md) 联网查询的结果取自 2026-10-29 当天。
>
> 排版提示：个别较宽的输出表是按 `options(width = 110)` 排版的。你的 R 控制台默认只有 80 个字符宽，这些表会折成两段，内容一样；想得到同样的排版，先运行一次 `options(width = 110)`。

---

## 先读这三条限定

1. **工具本身没有在这里运行。** LAVA、SUPERGNOVA、HDL、scDRS、scPagwas、CellOracle、scTenifoldKnk、NicheNet、Geneformer 这些工具，我没有安装和运行。每一章用的是"教学简化版"：在玩具数据上把原理演示一遍，不是工具本身。命令行和函数名依据文档写出，没有验证，以原文和官方文档为准。真正运行过的只有 `coloc`、`susieR`、`survival` 和联网查询。
2. **玩具数据是模拟的，真相已知。** 好处是能检验方法对不对；代价是里面的数字不能当作结论搬到 CHIP 研究里。
3. **报告里的 2025–2026 年文献结论、UK Biobank 字段与费用，我没有逐条核实。** 核实过的部分（GWAS 编号、PubMed 检索、Open Targets 注释）集中放在 [10 数据与注释](10-data-druggability-aac.md)，并列出了发现的几处与报告不一致的地方。

---

## 怎么读这套讲义

| 你想做什么 | 读哪几篇 |
|---|---|
| 先弄懂这份报告在做什么 | 本篇（[00 总览](00-overview-and-glossary.md)） |
| 想懂遗传学层（遗传相关、局部相关、共定位） | [02 遗传相关](02-genetic-correlation-ldsc-hdl.md)、[03 局部相关](03-local-genetic-correlation.md)、[04 共定位](04-colocalization-coloc-susie.md) |
| 想懂单细胞层 | [05 单细胞基础](05-single-cell-basics-pseudobulk.md)、[06 GWAS 映射到细胞](06-gwas-to-cells-magma-scdrs.md)、[07 程序评分](07-gene-program-scoring.md)；[08 虚拟扰动](08-in-silico-perturbation.md) 是进阶 |
| 想懂 UK Biobank 个体级分析 | [01 CHIP 与 VAF](01-chip-and-vaf.md)、[09 队列与共病](09-ukb-cohort-and-comorbidity.md) |
| 想查数据、可成药性、影像评分 | [10 数据与注释](10-data-druggability-aac.md) |
| 想把几步串起来、理解立项规则 | [11 流水线与立项规则](11-pipeline-and-go-nogo.md) |
| 想复习、考一考自己、查公式与文献 | [12 常见错误、自测题、速查卡与附录](12-pitfalls-quiz-and-reference.md) |
| 想学孟德尔随机化（MR） | [孟德尔随机化讲义](../mendelian-randomization.md) |

篇名里的记号：★ 是入门必读；☆ 是进阶内容，第一遍可以跳过。

---

## 先记住这几句话

1. 报告想解释一个"悖论"：骨丢失和血管钙化方向相反，却常常出现在同一批人身上。它的假说是：CHIP 突变的髓系细胞，在骨髓里和血管壁里开启了同一套"破骨样"程序，只是后果相反。
2. 要验证这个假说，需要三类证据：遗传学证据（两种表型是否共享遗传基础、是否共享同一个变异）、单细胞证据（风险是否富集在髓系细胞、两个组织的程序是否重叠）、队列证据（CHIP 携带者里骨和血管结局是否一起出现）。
3. 每一类证据都有一个常见的坑：遗传相关会被"正负抵消"掩盖；共定位会被"多个因果变异"骗；单细胞分析会把细胞当成独立样本；队列分析会被年龄混杂。本讲义每章先讲这个坑。
4. 纯计算只能做到"程序层面的比对"，要证明因果，报告自己也说小鼠实验不可省。本讲义不讲实验。

---

## 分篇目录

| 篇 | 主题 | 级别 |
|---|---|---|
| [00-overview-and-glossary.md](00-overview-and-glossary.md) | 总览、问题假说与术语速查 | ★ |
| [01-chip-and-vaf.md](01-chip-and-vaf.md) | CHIP 与 VAF：数据从哪来 | ★ |
| [02-genetic-correlation-ldsc-hdl.md](02-genetic-correlation-ldsc-hdl.md) | 遗传相关：LDSC 与 HDL | ★ |
| [03-local-genetic-correlation.md](03-local-genetic-correlation.md) | 局部遗传相关：LAVA 与 SUPERGNOVA | ★ |
| [04-colocalization-coloc-susie.md](04-colocalization-coloc-susie.md) | 共定位：coloc 与 SuSiE-coloc | ★ |
| [05-single-cell-basics-pseudobulk.md](05-single-cell-basics-pseudobulk.md) | 单细胞分析基础，以及"别把细胞当重复" | ★ |
| [06-gwas-to-cells-magma-scdrs.md](06-gwas-to-cells-magma-scdrs.md) | 把 GWAS 信号映射到细胞：MAGMA、scDRS 与 scPagwas | ★ |
| [07-gene-program-scoring.md](07-gene-program-scoring.md) | 基因程序评分与跨组织比对 | ★ |
| [08-in-silico-perturbation.md](08-in-silico-perturbation.md) | 虚拟扰动：CellOracle、scTenifoldKnk、NicheNet 与 Geneformer | ☆ |
| [09-ukb-cohort-and-comorbidity.md](09-ukb-cohort-and-comorbidity.md) | 个体级队列分析：CHIP × 骨密度 × 骨折，以及共病 | ★ |
| [10-data-druggability-aac.md](10-data-druggability-aac.md) | 数据与注释：GWAS Catalog、PubMed、Open Targets、影像 AAC 评分 | ☆ |
| [11-pipeline-and-go-nogo.md](11-pipeline-and-go-nogo.md) | 把它们串起来：三阶段流水线与"三取二"立项规则 | ★ |
| [12-pitfalls-quiz-and-reference.md](12-pitfalls-quiz-and-reference.md) | 常见错误、自测题、一页速查卡与附录 | ★ |

---

## 0. 术语速查（先扫一眼，看不懂再回来查）

| 术语 | 一句话解释 |
|---|---|
| CH、CHIP | 克隆性造血：造血干细胞带上体细胞突变并扩增成克隆。CHIP 指没有血液肿瘤的那一类，常用阈值是 VAF 不低于 2% |
| DNMT3A、TET2、ASXL1 | 最常见的三个 CHIP 驱动基因 |
| VAF | 变异等位基因比例：支持突变的读段数 ÷ 总读段数 |
| 大克隆、小克隆 | 报告里用 VAF 不低于 0.1 区分大克隆和小克隆 |
| BMD、eBMD、DXA | 骨密度；eBMD 是跟骨超声估计的骨密度；DXA 是双能 X 线测的骨密度 |
| CAC、AAC | 冠状动脉钙化、腹主动脉钙化 |
| 汇总统计 | GWAS 公开的结果表，每个 SNP 一行：beta、SE、p 等（详见 [MR 讲义](../mendelian-randomization.md) 第 5 章） |
| z 分数 | beta ÷ SE。z 的平方是 χ² 统计量 |
| LD、LD 分数 | 相邻 SNP 一起遗传；LD 分数是一个 SNP 与周围 SNP 的 r² 之和，衡量它"代表"了多少邻居 |
| 遗传度 h² | 性状的方差中，由 SNP 的加性效应解释的比例 |
| 遗传协方差、遗传相关 rg | 两个性状的遗传效应一起变化的程度；rg 是把它按两个遗传度标准化后的相关系数 |
| 全局、局部遗传相关 | 全局：整个基因组合起来；局部：某一小段基因组 |
| 有效样本量 | 病例对照研究换算成"相当于多少个均衡样本"：4 ÷ (1/病例 + 1/对照) |
| 共定位 | 判断两个性状在同一区域的信号，是否来自同一个因果变异 |
| PP.H0 到 PP.H4 | 共定位的五个假设的后验概率 |
| 可信集 | 一小撮 SNP，包含因果变异的概率不低于某个值（如 95%） |
| 计数矩阵 | 单细胞数据：行是细胞，列是基因，格子里是读到的分子数 |
| 伪重复 | 把同一个供体的很多细胞当成互相独立的重复，导致假阳性大量出现 |
| 伪批量 | 先把同一个供体的细胞汇总，再在供体层面做比较 |
| MAGMA | 把 SNP 水平的 GWAS 信号合并成基因水平分数的工具 |
| scDRS | 给每个单细胞打一个"疾病相关分数"的方法 |
| 基因程序、程序评分 | 一组一起升降的基因；程序评分是每个细胞里这组基因的"活跃程度" |
| GRN | 基因调控网络：谁调控谁 |
| 虚拟扰动 | 用计算模型预测基因被敲除或抑制后，细胞状态会怎样变化 |
| 线性基线 | 一个简单的线性预测，作为对照；复杂模型必须打赢它才算有价值 |
| 配体活性 | NicheNet 的输出：哪个配体最能解释受体细胞里观察到的基因变化 |
| 观察/期望比（O/E） | 共病分析：观察到的同时出现次数 ÷ 假设两件事独立时的期望次数 |
| Go/No-Go | 预实验结束时"继续还是放弃"的判定点 |

---

## 1. ★ 这份报告在做什么：问题、假说和技术地图

> 这一章带走什么：报告的核心问题和假说是什么；三阶段分析链各用了哪些技术；这些技术各自能回答什么、不能回答什么。

### 1.1 问题和假说

**观察。** 临床上，骨密度低的人，血管钙化往往反而更多，这叫"钙化悖论"。它很容易让人想到"钙从骨里流失，又沉积到血管里"，于是猜"骨丢失导致血管钙化"。

**一个负面结果。** 报告引用的文献说，用 MR 检验"骨密度 → 冠状动脉钙化"，两次都是阴性。也就是说，两者看起来不是因果关系。

**报告的假说。** 既然不是因果，那就是同一个上游的两个并行后果。衰老相关的克隆性造血（CHIP）使髓系细胞获得一套共同的"破骨样、促钙化"分泌程序：

- 在骨髓里，这套程序表现为破骨增强、骨丢失；
- 在血管壁里，这套程序表现为基质降解和代偿性钙化。

这个假说有一个好处：**可以证伪**。如果两个组织里的程序不重叠，假说就错了，而这本身也是一个可发表的结果。

### 1.2 三阶段分析链

下表是报告推荐的分析链（MR 的部分略去，另有讲义）。

| 阶段 | 步骤 | 技术 | 数据 | 本讲义 |
|---|---|---|---|---|
| 第一阶段（4–6 周，零成本） | 遗传相关 | LDSC、HDL | 公开 GWAS 汇总统计 | [02 遗传相关](02-genetic-correlation-ldsc-hdl.md) |
| | 局部遗传相关 | LAVA、SUPERGNOVA | 同上 | [03 局部相关](03-local-genetic-correlation.md) |
| | 共定位 | coloc、SuSiE-coloc | 同上 | [04 共定位](04-colocalization-coloc-susie.md) |
| | 风险细胞定位 | scDRS、scPagwas | 汇总统计 + 单细胞图谱 | [06 GWAS 映射到细胞](06-gwas-to-cells-magma-scdrs.md) |
| | 程序评分与比对 | 基因集评分 | CHIP 外周血、斑块、骨髓单细胞数据 | [07 程序评分](07-gene-program-scoring.md) |
| 第二阶段（3–6 个月） | 个体级分层分析、共病分析 | 回归、Cox 模型 | UK Biobank Tier 1 | [09 队列与共病](09-ukb-cohort-and-comorbidity.md) |
| | 可成药性注释 | Open Targets、drug2cell | 公共数据库 | [10 数据与注释](10-data-druggability-aac.md) |
| 第三阶段（可选） | 影像 AAC 评分 | 影像机器学习 | UK Biobank 侧位脊柱 DXA 影像 | [10 数据与注释](10-data-druggability-aac.md) |
| 贯穿 | 虚拟扰动 | CellOracle、scTenifoldKnk、NicheNet、Geneformer | 单细胞数据 | [08 虚拟扰动](08-in-silico-perturbation.md) |

### 1.3 每类证据能回答什么、不能回答什么

| 证据 | 能回答 | 不能回答 |
|---|---|---|
| 遗传相关、局部相关 | 两种表型是否共享遗传基础 | 共享的是不是同一个变异；是不是因果 |
| 共定位 | 某个区域里两个信号是否来自同一个变异 | 这个变异通过哪个基因、哪类细胞起作用 |
| 风险细胞定位（scDRS） | 遗传风险偏向哪类细胞 | 这类细胞里的哪个通路在起作用 |
| 程序评分与比对 | 两个组织里是否出现同一套基因程序 | 程序是原因还是结果 |
| 个体级队列分析 | CHIP 携带者里骨和血管结局是否一起出现 | 为什么一起出现 |
| 虚拟扰动 | 计算上抑制某个节点，程序会怎样变 | 真实生物学里会不会这样（要实验验证） |

报告自己指出了最关键的一个缺口：**没有任何公开数据集同时具备"人骨组织或骨髓 + CHIP 基因型"**。所以纯计算路径只能做到"程序层面的跨组织比对"，要做到因果，小鼠骨髓移植实验不可省。

### 1.4 立项规则：三取二

报告在第 6 周设了 Go/No-Go 判定点，以下三条满足两条则立项：

1. CHIP（或其任一亚型）与 eBMD、AAC、CAC 中至少两个表型存在显著的局部遗传相关或共定位信号；
2. scDRS 显示至少一个钙化相关 GWAS 的风险富集在髓系细胞上；
3. 破骨样程序评分在 CHIP 携带者单核细胞与钙化斑块巨噬细胞中同方向升高。

三条全部不满足，就不花 UK Biobank 的钱。[11 流水线与立项规则](11-pipeline-and-go-nogo.md) 会用几行 R 算一算这条规则的统计性质。

### 1.5 读报告时要留意的三处

这三处是我在核对时发现的，详见 [10 数据与注释](10-data-druggability-aac.md)：

1. 报告里的"HDL"指 **High-Definition Likelihood**（高清似然法），不是高密度脂蛋白胆固醇。
2. 报告把 CHIP 汇总统计（GCST90102618 到 622）的出处写成"Kessler 等，Nat Genet 2022"。GWAS Catalog 的元数据显示，这批数据来自 **Kar 等的 Nat Genet 2022**（PMID 35835912）。Kessler 等的那篇是 Nature 2022 的另一项工作。
3. 报告说"版本 B（ASXL1）遗传学层会缺一块"。但 2025 年那套多表型 CH GWAS（PMID 39948438，编号 GCST90435341 到 376）里**有 ASXL1**（GCST90435344、356、368）。

### 小结

- 报告的假说是"同一套破骨样程序在骨和血管里的并行后果"，它可以被证伪。
- 技术分三层：遗传学、单细胞、队列；每层有自己的坑。
- 纯计算做不到因果，只能做程序层面的比对。

---

[目录](README.md) · [下一篇：01 · CHIP 与 VAF：数据从哪来 →](01-chip-and-vaf.md)
