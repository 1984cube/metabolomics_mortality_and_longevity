[← 上一篇：11 · 把它们串起来：三阶段流水线与“三取二”立项规则](11-pipeline-and-go-nogo.md) · [目录](README.md)

---

# 12 · 常见错误、自测题、一页速查卡与附录

## 13. ★ 常见错误

| 技术 | 错误 | 后果 | 怎么避免 |
|---|---|---|---|
| LDSC、HDL | 忽略有效样本量和尺度 | 遗传度不可比、标准误大 | 用有效样本量；二分类性状换算尺度 |
| LDSC、HDL | 只看全局 rg | 正负抵消被漏掉 | 加局部方法 |
| LDSC、HDL | 用不匹配的人群的 LD 分数 | 偏倚 | 人群一致 |
| 局部 rg | 不检验局部遗传度就估计 rg | rg 不稳定 | 先过滤 |
| 局部 rg | 不校正区块数 | 假阳性 | Bonferroni 或 FDR |
| coloc | 等位基因没对齐、基因组版本不一致 | 结果无意义 | 先协调 |
| coloc | 多信号区域仍用单变量 coloc | 误判 H3（[04 共定位](04-colocalization-coloc-susie.md) 的例子里 PP.H3 = 1） | 用 SuSiE-coloc |
| coloc | SuSiE 的 LD 矩阵与人群不匹配 | 假的可信集 | 同人群 LD |
| coloc | 把低 PP.H4 当成"没有共定位" | 误读功效不足 | 同时看 PP.H0、PP.H1、PP.H2 |
| 单细胞 | 把细胞当重复 | 假阳性（[05 单细胞基础](05-single-cell-basics-pseudobulk.md) 的例子里 67%） | 供体层面分析 |
| 单细胞 | 直接比较不同数据集的原始分数 | 批次、平台造成的差异 | 比较方向和排序 |
| scDRS | 用朴素的随机对照 | 所有细胞类型都"显著" | 表达量匹配的对照 |
| scDRS | 图谱里没有目标细胞类型 | 漏掉真正的靶细胞 | 先确认图谱覆盖 |
| 程序比对 | 只报 Spearman，用简单置换 | 共同基线造成假相关 | 前 N 个基因重叠 + 分箱置换 |
| 程序评分 | 评分基因包含被扰动的候选节点 | 循环论证 | 评分基因集排除候选节点 |
| 虚拟扰动 | 不与线性基线并排 | 高估复杂模型 | 同时报告基线 |
| 虚拟扰动 | 没有方向先验的网络 | 把方向弄反 | 底图网络 |
| 队列 | 不调整年龄 | 效应放大近 6 倍 | 调整年龄和性别等 |
| 队列 | 一个亚组显著、一个不显著就说有差别 | 误判 | 交互项或对比检验 |
| 共病 | O/E > 1 就当成共享程序的证据 | 共同原因也能造成 | 与非携带者比较、看剂量反应 |
| Go/No-Go | 事后定阈值 | 假阳性 | 事先写进方案 |
| 数据 | 把出处写错 | 引用错误 | 用 API 核对编号和 PMID |

---

## 14. ★ 自测题

先自己答，再看后面的答案。

1. VAF = 0.10 的杂合突变，大约有多少比例的细胞带这个突变？
2. 为什么外显子测序对小克隆不敏感？
3. LDSC 的截距和斜率分别代表什么？
4. 为什么 LDSC 估计的遗传相关对样本重叠稳健？
5. 病例 2,041、对照 173,918，有效样本量大约是多少？这意味着什么？
6. 全局 rg 约等于 0，能不能说两个性状没有共享？
7. 局部遗传相关方法为什么要先检验局部遗传度？
8. coloc 里 PP.H3 很高，一定是"两个性状因果变异不同"吗？
9. SuSiE-coloc 比 coloc 多做了什么？
10. 比较携带者和对照的单细胞表达，统计的重复单位是什么？为什么？
11. scDRS 为什么要用匹配的对照基因集？
12. 跨组织比对时，为什么简单置换会给出假阳性？
13. 虚拟扰动为什么要并排报告线性基线？
14. 共病分析的 O/E 比为什么要在年龄 × 性别的分层里算，并和非携带者比？
15. "三取二"规则为什么要事先定好每条判据的阈值？

### 参考答案

1. 约 20%。带突变的细胞比例约等于 2 × VAF。
2. 小克隆的 VAF 低，在有限的测序深度下，支持突变的读段数太少，达不到检出阈值。所以数据里的小克隆是被"截尾"的，"非携带者"里混着看不见的小克隆。
3. 斜率 = N h² ÷ M，给出遗传度；截距 = 1 + N a，反映人群分层等混杂，没有混杂时约等于 1。
4. 样本重叠只改变跨性状回归的截距（ρ N_s ÷ √(N1 N2)），不改变斜率；rg 由斜率和两个遗传度算出。
5. 4 ÷ (1/2,041 + 1/173,918) ≈ 8,069。相当于一个只有约 8,000 人的均衡样本，LDSC 的功效会很低，遗传度和遗传相关的估计都不稳。
6. 不能。正负相关的区域可以互相抵消，使全局 rg 接近 0。要用局部方法逐区块检验。
7. 没有局部遗传信号时，局部 rg 的分母接近 0，估计极不稳定。先检验局部遗传度，只在两个性状都有信号的区块里估计 rg。
8. 不一定。区域里有多个独立信号、而真正共享的是较弱的次要信号时，单一因果变异的 coloc 也会给出高的 PP.H3（[04 共定位](04-colocalization-coloc-susie.md) 5.6 节的例子）。要用 SuSiE-coloc 按可信集配对。
9. 先用 SuSiE 把每个性状的信号拆成若干可信集，再对每一对可信集做共定位，从而允许一个区域里有多个因果变异。需要匹配的 LD 矩阵。
10. 供体（人）。同一个供体的细胞共享这个供体自己的基线，不是独立样本。把细胞当重复，假阳性会大量出现（[05 单细胞基础](05-single-cell-basics-pseudobulk.md) 6.5 节的例子里 67%）。
11. 疾病基因集往往表达量偏高。与随机基因比，所有细胞类型的分数都会偏高，没有分辨力。匹配表达量的对照才能反映"这个基因集在某类细胞里是否额外偏高"。
12. 两个组织共用同一批基因，基线表达量带来的、被标准化放大的系统性偏移在两边是共同的，造出假的相关。简单置换把基因随机打乱，破坏了这种结构，所以看不出来。应在同一表达量分档内置换。
13. 线性基线经常已经够好（如根转录因子的情形），复杂模型只有明显优于基线才有价值；没有基线，就无法判断复杂模型有没有真正学到东西。
14. 年龄会让两件事天然结伴出现，分层把它拿掉；吸烟、肾功能等共同原因在非携带者里也会让 O/E 偏离 1，所以要和非携带者比，看携带者是否额外结伴。
15. 事后定阈值会让假阳性不可控；而且判据 1 里有大量隐含的多重检验，不事先写清楚校正范围，立项判断就失去统计意义。

---

## 15. ★ 一页速查卡

| 技术 | 一句话 | 先看哪个数 | 最常见的坑 |
|---|---|---|---|
| VAF | 带突变的细胞比例约 2 × VAF | VAF 与测序深度 | 小克隆被漏掉 |
| LDSC | χ² 对 LD 分数回归：斜率给遗传度，截距给混杂 | 截距、h² 的 z 分数 | 有效样本量太小 |
| HDL | LDSC 的全似然版本，方差更小 | 与 LDSC 方向是否一致 | 需要专门的 LD 参考面板 |
| LAVA、SUPERGNOVA | 逐区块估计局部遗传相关 | 局部遗传度是否显著、区块数 | 多重检验 |
| coloc | 两个信号是否来自同一个变异 | PP.H4（不低于 0.8）、PP.H3 | 多信号区域 |
| SuSiE-coloc | 按可信集配对做共定位 | 每对可信集的 PP.H4 | LD 不匹配 |
| 单细胞差异分析 | 重复单位是供体 | 每组供体数、最小可检出效应 | 把细胞当重复 |
| scDRS | 匹配对照下的细胞疾病分数 | 细胞类型层面的 z 与显著比例 | 朴素对照 |
| 程序评分 | 基于排名的基因集评分 | 供体层面的 p 值 | 评分基因含候选节点 |
| 跨组织比对 | 相关 + 前 N 重叠 + 分箱置换 | 分箱置换 p | 共同基线的假相关 |
| 虚拟扰动 | 在计算机里敲除，看评分怎么变 | 与线性基线的差值 | 无先验、不比基线 |
| 个体级分析 | 调整年龄后按基因和克隆大小分层 | 亚组样本量、最小可检出差异 | 年龄混杂 |
| 共病 O/E | 观察/期望比，在分层里算，与非携带者比 | 携带者 O/E 与非携带者 O/E | 共同原因 |
| 三取二 | 三条判据满足两条就立项 | 每条的假阳性率和功效 | 事后定阈值 |

---

## 附录 A　公式一览

```text
VAF 与细胞比例        带突变的细胞比例 = 2 × VAF（常染色体杂合突变）
有效样本量            N_eff = 4 ÷ (1/N_病例 + 1/N_对照)
LD 分数               ℓ_j = Σ_k r²_jk（校正：r² − (1 − r²)/(N − 2)）
LDSC 单性状           E[χ²_j] = N h² ℓ_j / M + N a + 1
LDSC 跨性状           E[z1_j z2_j] = √(N1 N2) ρ_g ℓ_j / M + ρ N_s / √(N1 N2)
遗传相关              rg = ρ_g ÷ √(h²₁ h²₂)
近似贝叶斯因子        r = W ÷ (W + V)，z = beta ÷ SE，log ABF = 0.5 × (log(1 − r) + r z²)，W = (0.15 sdY)²
coloc 五个假设        H0 = 1；H1 = p1 Σ ABF1；H2 = p2 Σ ABF2；
                      H3 = p1 p2 (Σ ABF1 × Σ ABF2 − Σ ABF1 ABF2)；H4 = p12 Σ ABF1 ABF2；再归一化
UCell                 U = 1 − (Σ r_i − n(n + 1)/2) ÷ (n × maxRank)
scDRS（教学版）       z_c = (s_c − 对照均值) ÷ 对照标准差；p_c = (1 + #{对照 ≥ s_c}) ÷ (B + 1)
两组供体的功效        非中心参数 = d √(n1 n0 ÷ (n1 + n0))，自由度 = n1 + n0 − 2
观察/期望比           O/E = Σ 观察的同时发生数 ÷ Σ（分层人数 × P(事件A) × P(事件B)）
三取二                P(至少两条满足) = p1 p2 (1 − p3) + p1 (1 − p2) p3 + (1 − p1) p2 p3 + p1 p2 p3
```

---

## 附录 B　符号表

| 符号 | 含义 |
|---|---|
| VAF | 变异等位基因比例 |
| N、M | 样本量、SNP 个数 |
| h² | 遗传度 |
| ρ_g | 遗传协方差 |
| rg | 遗传相关 |
| ℓ_j | 第 j 个 SNP 的 LD 分数 |
| a | LDSC 里衡量混杂的参数 |
| N_s | 两份 GWAS 重叠的样本数 |
| z、χ² | z 分数、卡方统计量（z²） |
| ABF | 近似贝叶斯因子 |
| W、V | 先验效应方差、效应估计的方差 |
| p1、p2、p12 | coloc 的先验概率 |
| PP.H0 到 PP.H4 | coloc 五个假设的后验概率 |
| sdY | 性状的标准差 |
| U | UCell 评分 |
| O/E | 观察/期望比 |
| HR | 风险比 |

---

## 附录 C　文献

下面的题名是简写，检索时请用作者、期刊和年份。

CHIP 与数据：

- Jaiswal S 等. Age-related clonal hematopoiesis associated with adverse outcomes. N Engl J Med, 2014.
- Genovese G 等. Clonal hematopoiesis and blood-cancer risk inferred from blood DNA sequence. N Engl J Med, 2014.
- Steensma DP 等. Clonal hematopoiesis of indeterminate potential and its distinction from myelodysplastic syndromes. Blood, 2015.
- Vlasschaert C 等. A practical approach to curate clonal hematopoiesis of indeterminate potential in human genetic data sets. Blood, 2023.
- Kar SP 等. Genome-wide analyses of 200,453 individuals yield new insights into the causes and consequences of clonal hematopoiesis. Nat Genet, 2022（PMID 35835912）.
- Kessler MD 等. Common and rare variant associations with clonal haematopoiesis phenotypes. Nature, 2022（PMID 36450978）.
- Wen S 等. Comparative analysis of the Mexico City Prospective Study and the UK Biobank identifies ancestry-specific effects on clonal hematopoiesis. Nat Genet, 2025（PMID 39948438）.
- Morris JA 等. An atlas of genetic influences on osteoporosis in humans and mice. Nat Genet, 2019.
- Kavousi M 等. Multi-ancestry genome-wide study identifies effector genes and druggable pathways for coronary artery calcification. Nat Genet, 2023.
- Sethi A 等. Calcification of the abdominal aorta is an under-appreciated cardiovascular disease risk factor in the general population. Front Cardiovasc Med, 2022.

遗传相关与局部相关：

- Bulik-Sullivan BK 等. LD Score regression distinguishes confounding from polygenicity in genome-wide association studies. Nat Genet, 2015.
- Bulik-Sullivan B 等. An atlas of genetic correlations across human diseases and traits. Nat Genet, 2015.
- Ning Z, Pawitan Y, Shen X. High-definition likelihood inference of genetic correlations across human complex traits. Nat Genet, 2020.
- Werme J 等. An integrated framework for local genetic correlation analysis. Nat Genet, 2022.
- Zhang Y 等. SUPERGNOVA: local genetic correlation analysis reveals heterogeneous etiologic sharing of complex traits. Genome Biol, 2021.

共定位：

- Giambartolomei C 等. Bayesian test for colocalisation between pairs of genetic association studies using summary statistics. PLoS Genet, 2014.
- Wakefield J. Bayes factors for genome-wide association studies: comparison with P-values. Genet Epidemiol, 2009.
- Wallace C. A more accurate method for colocalisation analysis allowing for multiple causal variants. PLoS Genet, 2021.
- Wang G 等. A simple new approach to variable selection in regression, with application to genetic fine mapping. J R Stat Soc B, 2020.
- Foley CN 等. A fast and efficient colocalization algorithm for identifying shared genetic risk factors across multiple traits. Nat Commun, 2021.

单细胞：

- Wolf FA 等. SCANPY: large-scale single-cell gene expression data analysis. Genome Biol, 2018.
- Hao Y 等. Integrated analysis of multimodal single-cell data. Cell, 2021.
- Squair JW 等. Confronting false discoveries in single-cell differential expression. Nat Commun, 2021.
- Zimmerman KD 等. A practical solution to pseudoreplication bias in single-cell studies. Nat Commun, 2021.
- de Leeuw CA 等. MAGMA: generalized gene-set analysis of GWAS data. PLoS Comput Biol, 2015.
- Zhang MJ 等. Polygenic enrichment distinguishes disease associations of individual cells in single-cell RNA-seq data. Nat Genet, 2022.
- Ma Y 等. Polygenic regression uncovers trait-relevant cellular contexts through pathway activation transformation of single-cell RNA sequencing data. Cell Genom, 2023.
- Andreatta M, Carmona SJ. UCell: robust and scalable single-cell gene signature scoring. Comput Struct Biotechnol J, 2021.
- Aibar S 等. SCENIC: single-cell regulatory network inference and clustering. Nat Methods, 2017.
- Tirosh I 等. Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq. Science, 2016.

虚拟扰动：

- Kamimoto K 等. Dissecting cell identity via network inference and in silico gene perturbation. Nature, 2023.
- Osorio D 等. scTenifoldKnk: an efficient virtual knockout tool for gene function predictions via single-cell gene regulatory network perturbation. Patterns, 2022.
- Browaeys R, Saelens W, Saeys Y. NicheNet: modeling intercellular communication by linking ligands to target genes. Nat Methods, 2020.
- Theodoris CV 等. Transfer learning enables predictions in network biology. Nature, 2023.
- Ahlmann-Eltze C, Huber W, Anders S. Deep learning-based predictions of gene perturbation effects do not yet outperform simple linear baselines. Nat Methods, 2025.

队列、注释与影像：

- Gelman A, Stern H. The difference between "significant" and "not significant" is not itself statistically significant. Am Stat, 2006.
- Ochoa D 等. The next-generation Open Targets Platform: reimagined, redesigned, rebuilt. Nucleic Acids Res, 2023.
- Kauppila LI 等. New indices to classify location, severity and progression of calcific lesions in the abdominal aorta: a 25-year follow-up study. Atherosclerosis, 1997.

---

## 核对记录与版本说明

**我核对过的（2026-10-29，联网）：**

1. 上传的两份 PDF 逐字节相同（MD5 一致），只有一份报告。
2. 报告里 10 个 GWAS 编号都用 GWAS Catalog 的 REST API 查过：性状和样本量与报告一致。
3. 编号对应的论文用 PubMed 查过：GCST90102618 到 622 是 Kar 等 Nat Genet 2022（PMID 35835912），不是报告写的 Kessler 等；Kessler 等的是 Nature 2022（PMID 36450978）。
4. ASXL1 的汇总统计在 GCST90435344、356、368 中存在（PMID 39948438，Nat Genet 2025）。
5. CHIP 汇总统计的真实文件头已取回，确认是 BOLT-LMM 的结果，坐标为 GRCh37。
6. [10 数据与注释](10-data-druggability-aac.md) 的 PubMed 检索、Open Targets 查询都真实运行过，结果就是表里显示的。
7. 所有 R 代码块都从本讲义的源文件自动嵌入，在 R 4.3.3 里按章节运行；`coloc` 5.2.3、`susieR` 0.14.2、`jsonlite` 2.0.0。我手写的 coloc 实现与 `coloc` 包的结果完全一致。

**我没有核对的：**

- 报告里 2025–2026 年发表的机制论文的具体结论（JCI 2026、Bone 2026、J Immunol 2026、Nat Methods 2026 的基准研究等）；
- UK Biobank 的字段编号、费用、档位和审批周期；
- 报告对竞争格局的判断（只核对了一组检索，见 [10 数据与注释](10-data-druggability-aac.md)）；
- LAVA、SUPERGNOVA、HDL、scDRS、scPagwas、CellOracle、scTenifoldKnk、NicheNet、Geneformer、drug2cell 这些工具本身。本讲义对它们的描述来自原文的公开信息，没有运行；
- 讲义里 LDSC 的命令行示例没有运行。

**可复现性：** 涉及随机数的结果（所有模拟）在不同的 R 版本或不同的随机种子下，小数位可能略有差别，结论不变。[02 遗传相关](02-genetic-correlation-ldsc-hdl.md) 的 LDSC 模拟单次估计噪声较大，所以正文里用了重复和标准误来呈现。

---

[← 上一篇：11 · 把它们串起来：三阶段流水线与“三取二”立项规则](11-pipeline-and-go-nogo.md) · [目录](README.md)
