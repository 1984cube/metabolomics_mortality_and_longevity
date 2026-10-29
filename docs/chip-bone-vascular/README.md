# CHIP × 骨–血管钙化课题：生信技术讲义（分篇目录）

这是《CHIP × 骨–血管钙化悖论：课题可行性核实报告》的配套生信讲义（孟德尔随机化除外；MR 请看 [`../mendelian-randomization.md`](../mendelian-randomization.md)）。

为了方便按专题阅读，讲义拆成了 13 篇独立文件。每一篇都可以单独打开阅读，里面的 R 代码也可以单独运行（代码在 R 4.3.3 中实际跑过，输出自动嵌入）。

## 分篇列表

| 序号 | 文件 | 主题 | 行数 |
|---:|---|---|---:|
| 00 | [`00-overview-and-glossary.md`](00-overview-and-glossary.md) | 00 · 总览、问题假说与术语速查 | 176 |
| 01 | [`01-chip-and-vaf.md`](01-chip-and-vaf.md) | 01 · CHIP 与 VAF：数据从哪来 | 131 |
| 02 | [`02-genetic-correlation-ldsc-hdl.md`](02-genetic-correlation-ldsc-hdl.md) | 02 · 遗传相关：LDSC 与 HDL | 312 |
| 03 | [`03-local-genetic-correlation.md`](03-local-genetic-correlation.md) | 03 · 局部遗传相关：LAVA 与 SUPERGNOVA | 237 |
| 04 | [`04-colocalization-coloc-susie.md`](04-colocalization-coloc-susie.md) | 04 · 共定位：coloc 与 SuSiE-coloc | 248 |
| 05 | [`05-single-cell-basics-pseudobulk.md`](05-single-cell-basics-pseudobulk.md) | 05 · 单细胞分析基础，以及“别把细胞当重复” | 184 |
| 06 | [`06-gwas-to-cells-magma-scdrs.md`](06-gwas-to-cells-magma-scdrs.md) | 06 · 把 GWAS 信号映射到细胞：MAGMA、scDRS 与 scPagwas | 147 |
| 07 | [`07-gene-program-scoring.md`](07-gene-program-scoring.md) | 07 · 基因程序评分与跨组织比对 | 171 |
| 08 | [`08-in-silico-perturbation.md`](08-in-silico-perturbation.md) | 08 · 虚拟扰动：CellOracle、scTenifoldKnk、NicheNet 与 Geneformer | 222 |
| 09 | [`09-ukb-cohort-and-comorbidity.md`](09-ukb-cohort-and-comorbidity.md) | 09 · 个体级队列分析：CHIP × 骨密度 × 骨折，以及共病 | 253 |
| 10 | [`10-data-druggability-aac.md`](10-data-druggability-aac.md) | 10 · 数据与注释：GWAS Catalog、PubMed、Open Targets、影像 AAC 评分 | 207 |
| 11 | [`11-pipeline-and-go-nogo.md`](11-pipeline-and-go-nogo.md) | 11 · 把它们串起来：三阶段流水线与“三取二”立项规则 | 153 |
| 12 | [`12-pitfalls-quiz-and-reference.md`](12-pitfalls-quiz-and-reference.md) | 12 · 常见错误、自测题、一页速查卡与附录 | 228 |

## 按目的选读

- **只想先看懂报告在干什么**：读 [`00-overview-and-glossary.md`](00-overview-and-glossary.md)
- **做第一阶段遗传学预实验**：读 `02` → `03` → `04`
- **做第一阶段单细胞预实验**：读 `05` → `06` → `07`（进阶再看 `08`）
- **准备 UK Biobank 个体级分析**：读 `01` → `09`
- **查公开数据、靶点可成药性、影像评分**：读 `10`
- **看整体流水线与第 6 周 Go/No-Go 判据**：读 `11`
- **考一考自己、查常见错误与公式速查卡**：读 `12`
