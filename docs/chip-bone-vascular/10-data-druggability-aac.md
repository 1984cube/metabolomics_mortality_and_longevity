[← 上一篇：09 · 个体级队列分析：CHIP × 骨密度 × 骨折，以及共病](09-ukb-cohort-and-comorbidity.md) · [目录](README.md) · [下一篇：11 · 把它们串起来：三阶段流水线与“三取二”立项规则 →](11-pipeline-and-go-nogo.md)

---

# 10 · 数据与注释：GWAS Catalog、PubMed、Open Targets、影像 AAC 评分

> 这一章带走什么：怎么查报告里的 GWAS 编号、怎么下载和读汇总统计文件；怎么用检索量化"没人做过"这个说法；怎么给候选基因做可成药性注释；腹主动脉钙化的影像评分怎么评价。这一章的查询是 2026-10-29 当天联网运行的。

### 11.1 查 GWAS Catalog：编号与元数据

GWAS Catalog 的 REST API（v2）可以按编号查研究的元数据。下面把报告里的编号逐个查一遍：

```r
library(jsonlite)
get_study <- function(id) {                           # GWAS Catalog REST API v2：按编号查研究的元数据
  d <- fromJSON(paste0("https://www.ebi.ac.uk/gwas/rest/api/v2/studies/", id))
  data.frame(编号 = id, 性状 = d$disease_trait, 样本 = d$initial_sample_size, PMID = d$pubmed_id)
}
ids <- c("GCST90102618", "GCST90102619", "GCST90102620", "GCST90102621", "GCST90102622",
         "GCST006979", "GCST90029004", "GCST90278456", "GCST90134614",
         "GCST90435344", "GCST90435356", "GCST90435368")
tab <- do.call(rbind, lapply(ids, get_study))
tab$样本 <- sub(" individuals", "", gsub(" ancestry", "", tab$样本))
tab$样本 <- ifelse(nchar(tab$样本) > 46, paste0(substr(tab$样本, 1, 43), "..."), tab$样本)
print(tab, row.names = FALSE)
```

```text
         编号                                   性状                                           样本     PMID
 GCST90102618         Clonal hematopoiesis (overall) 10,203 European cases, 173,918 European con... 35835912
 GCST90102619 Clonal hematopoiesis (DNMT3A mutation) 5,185 European cases, 173,918 European cont... 35835912
 GCST90102620   Clonal hematopoiesis (TET2 mutation) 2,041 European cases, 173,918 European cont... 35835912
 GCST90102621     Clonal hematopoiesis (large clone) 4,049 European cases, 173,918 European cont... 35835912
 GCST90102622     Clonal hematopoiesis (small clone) 6,154 European cases, 173,918 European cont... 35835912
   GCST006979              Heel bone mineral density                                426,824 British 30598549
 GCST90029004              Heel bone mineral density                               583,314 European 29892013
 GCST90278456          Coronary artery calcification                                28,655 European 37770635
 GCST90134614  Abdominal aortic calcification levels                           31,786 White British 36277789
 GCST90435344  Clonal hematopoiesis (ASXL1 mutation)             132,546 Hispanic or Latin American 39948438
 GCST90435356  Clonal hematopoiesis (ASXL1 mutation)                               388,769 European 39948438
 GCST90435368  Clonal hematopoiesis (ASXL1 mutation) 132,546 Hispanic or Latin American, 388,769... 39948438
```

**核对结果：**

1. **报告里的 GWAS 编号都能查到，性状、样本量与报告的描述一致。** CHIP 总体 10,203 例病例、173,918 例对照，DNMT3A-CH 5,185 例，TET2-CH 2,041 例，大克隆 4,049 例，小克隆 6,154 例；跟骨 BMD 426,824 人和 583,314 人；CAC 28,655 人；AAC 31,786 人，都对得上。
2. **出处与报告不一致。** 报告把 GCST90102618 到 622 写成"Kessler 等，Nat Genet 2022"。GWAS Catalog 显示它们对应 PMID 35835912，用 PubMed 查是 **Kar 等，Nat Genet 2022**（"Genome-wide analyses of 200,453 individuals yield new insights into the causes and consequences of clonal hematopoiesis"）。Kessler 等的那篇题为"Common and rare variant associations with clonal haematopoiesis phenotypes"，发表在 **Nature 2022（PMID 36450978）**，不是同一项研究。引用时要改。
3. **ASXL1 的汇总统计是有的。** 报告说"Kessler 2022 的公开汇总统计只有 DNMT3A-CH 和 TET2-CH，没有 ASXL1-CH，遗传学层会缺一块"。但 2025 年那套多表型 CH GWAS（PMID 39948438，Wen 等，Nat Genet 2025）在 GWAS Catalog 里有 ASXL1：欧裔 GCST90435356（388,769 人）、拉美裔 GCST90435344（132,546 人）、两者合并 GCST90435368。同一套数据里还有 PPM1D、TP53、SF3B1、SRSF2、GNB1、IDH2、JAK2 等。样本量一栏是参与分析的总人数，不是 ASXL1 病例数，病例数需要看原文。所以版本 B（ASXL1）在遗传学层上不缺数据。

### 11.2 下载和读汇总统计文件

每个编号在 FTP 上有一个目录。用 `curl` 只取压缩文件的开头，再解压看前几行，不必下载整个 392 MB：

```bash
curl -s -r 0-8000 \
  "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90102001-GCST90103000/GCST90102618/GCST90102618_buildGRCh37.tsv.gz" \
  | gunzip -c 2>/dev/null | head -3
```

输出（制表符分隔，下表按列整理）：

| variant_id | chromosome | base_pair_location | effect_allele | other_allele | effect_allele_frequency | INFO | beta | standard_error | CHISQ_BOLT_LMM | p_value |
|---|---|---|---|---|---|---|---|---|---|---|
| rs367896724 | 1 | 10177 | A | AC | 0.60091 | 0.467935 | 2.754e-4 | 0.00111752 | 0.0496341 | 0.82 |
| rs201106462 | 1 | 10352 | T | TA | 0.606925 | 0.447895 | -5.32636e-4 | 0.00115018 | 0.20858 | 0.65 |

三点提醒：

1. **beta 的尺度。** 有一列 `CHISQ_BOLT_LMM`，说明是 BOLT-LMM 线性混合模型的结果，对二分类的 CHIP，beta 是线性尺度（每多一个等位基因，患病概率的变化），不是 log OR。LDSC 的遗传度要注意"观察尺度 vs liability 尺度"，coloc 要按连续性状处理（见 [02 遗传相关](02-genetic-correlation-ldsc-hdl.md) 与 [04 共定位](04-colocalization-coloc-susie.md)）。
2. **基因组版本。** 文件名里的 `buildGRCh37` 说明坐标是 GRCh37。几套数据合并前，要确认版本一致，不一致要先转换（如 liftOver）。
3. **等位基因。** 做遗传相关、共定位之前，要对齐效应等位基因（见 [MR 讲义](../mendelian-randomization.md) 第 5 章）。

### 11.3 用检索量化"没人做过"

报告说："用 PubMed 同时检索 CHIP、骨与钙化，49 条结果里没有任何一篇把骨和血管放在同一项研究里。"这类说法可以用 NCBI E-utilities 复现：

```r
pubmed_count <- function(term) {                       # NCBI E-utilities：只取命中条数
  url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmode=json&retmax=0&term=",
                URLencode(term, reserved = TRUE))
  Sys.sleep(0.5)                                       # 不带 API key 时，每秒最多 3 次请求
  as.integer(fromJSON(url)$esearchresult$count)
}
ch <- '"clonal hematopoiesis"[tiab]'
q  <- c(`CH + 骨质疏松/骨密度/骨丢失` = paste(ch, 'AND (osteoporosis[tiab] OR "bone mineral density"[tiab] OR "bone loss"[tiab])'),
        `CH + 钙化/主动脉瓣狭窄`      = paste(ch, 'AND (calcification[tiab] OR "aortic stenosis"[tiab])'),
        `CH + 骨 + 钙化`              = paste(ch, 'AND (bone[tiab] OR osteoporosis[tiab]) AND calcification[tiab]'))
data.frame(检索 = names(q), 命中条数 = sapply(q, pubmed_count), row.names = NULL)
```

```text
                         检索 命中条数
1 CH + 骨质疏松/骨密度/骨丢失       12
2      CH + 钙化/主动脉瓣狭窄       21
3              CH + 骨 + 钙化        2
```

我的三条检索式命中 12、21 和 2 条（和报告的 49 条不同，因为检索式不同，数字本来就会变）。第三条命中的 2 篇是：JCI 2026 年关于"CHIP 激活巨噬细胞促钙化通路、促进主动脉瓣狭窄"的论文，以及 NEJM 2017 年"克隆性造血与动脉粥样硬化性心血管病风险"的论文，**都不是骨与血管的联合研究**，与报告的判断一致。

使用要点：命中条数只是起点，**必须读标题和摘要**才能判断；要记录检索式和检索日期；这类检索只能说"我没找到"，不能证明"没人做过"（预印本、会议摘要、未发表的工作都检索不到）。

### 11.4 可成药性注释：Open Targets 与 drug2cell

候选节点找出来之后，要问：**能不能被药物靶向？** Open Targets 平台汇总了靶点与疾病的关联证据，以及靶点的"可成药性"（tractability）：小分子、抗体等不同形式各处于什么阶段。下面用它的 GraphQL 接口查报告里的几个基因：

```r
ot <- function(query) fromJSON(system2("curl", c("-s", "-H", shQuote("Content-Type: application/json"),
        "-d", shQuote(toJSON(list(query = query), auto_unbox = TRUE)), "https://api.platform.opentargets.org/api/v4/graphql"),
        stdout = TRUE))$data
genes <- c("MMP9", "CTSK", "ACP5", "OSM", "IL20", "NFATC1", "TET2")
order_lv <- c("Approved Drug", "Advanced Clinical", "Phase 1 Clinical", "Structure with Ligand", "High-Quality Ligand",
              "High-Quality Pocket", "Druggable Family")                          # 从“最成熟”到“最早期”
best <- function(tr, mod) {
  lab <- tr$label[tr$modality == mod & tr$value]
  hit <- order_lv[order_lv %in% lab]
  if (length(hit)) hit[1] else "（无）"
}
rows <- lapply(genes, function(g) {
  id <- ot(sprintf('{ search(queryString: "%s", entityNames: ["target"], page: {index: 0, size: 1}) { hits { id } } }', g))$search$hits$id
  tr <- ot(sprintf('{ target(ensemblId: "%s") { tractability { label modality value } } }', id))$target$tractability
  data.frame(基因 = g, Ensembl = id, 小分子最成熟证据 = best(tr, "SM"), 抗体最成熟证据 = best(tr, "AB"))
})
do.call(rbind, rows)
```

```text
    基因         Ensembl      小分子最成熟证据    抗体最成熟证据
1   MMP9 ENSG00000100985     Advanced Clinical Advanced Clinical
2   CTSK ENSG00000143387     Advanced Clinical            （无）
3   ACP5 ENSG00000102575                （无）            （无）
4    OSM ENSG00000099985                （无）            （无）
5   IL20 ENSG00000162891                （无） Advanced Clinical
6 NFATC1 ENSG00000131196   High-Quality Ligand            （无）
7   TET2 ENSG00000168769 Structure with Ligand            （无）
```

读法：

- MMP9：小分子和抗体都有"Advanced Clinical"（进入后期临床）的证据；
- CTSK：小分子有后期临床证据；
- IL20：抗体有后期临床证据；
- ACP5、OSM：没有可成药性的证据；
- NFATC1、TET2：只有早期线索（"高质量配体"、"带配体的结构"）。

**这些标签只表示"可成药性证据的成熟程度"，不表示"有效"或"安全"。** 一个基因在后期临床里有药，不代表对骨和血管同时有益，更要看同时抑制两侧的假说。

报告还提到 drug2cell：把药物与靶点的关系映射到单细胞表达数据上，看一种药的靶点在哪些细胞里表达，用来判断药物可能作用于哪类细胞。本讲义没有运行这个工具，以原文为准。

### 11.5 影像 AAC 评分

腹主动脉钙化（AAC）是报告里"唯一的短板"：UK Biobank 里没有现成的 AAC 评分，要在 5 万张侧位脊柱 DXA 影像上自建。

**评分标准。** 常用的是 Kauppila 的 AAC-24：看腰椎 L1 到 L4 四个节段的腹主动脉前壁和后壁，每个位置按钙化长度打 0 到 3 分（0 无；1 小于 1/3；2 在 1/3 到 2/3 之间；3 大于 2/3），共 8 项，满分 24 分。

**自动评分模型怎么评价。** 模型从影像预测 AAC-24，要和人工评分比较，常看三个指标：

- Spearman 相关：排序是否一致；
- 二次加权 kappa：考虑"差几分"的一致性；
- 临床分组的一致率：如 0 分、1 到 5 分、6 分以上（示例分组）。

下面用模拟数据演示这三个指标怎么算：

```r
# Kauppila 的 AAC-24：腰椎 L1–L4 四个节段，每段的前壁、后壁各打 0–3 分（钙化长度：0 无；1 小于 1/3；2 在 1/3 到 2/3 之间；3 大于 2/3），共 8 项，满分 24
aac24 <- function(scores) { stopifnot(length(scores) == 8, all(scores %in% 0:3)); sum(scores) }
aac24(c(0, 0, 1, 0, 2, 1, 3, 2))

# 评估一个自动评分模型：和人工评分比较
set.seed(81)
n <- 500
manual <- pmin(24, rpois(n, 3) * (runif(n) < 0.6))                 # 人工评分（模拟）：很多人是 0 分
model  <- pmin(24, pmax(0, round(manual + rnorm(n, 0, 1.5))))      # 模型评分 = 人工评分 + 误差（模拟）
wkappa <- function(a, b, K = 25) {                                 # 二次加权 kappa（0–24 分共 25 档）
  W <- outer(0:(K - 1), 0:(K - 1), function(i, j) (i - j)^2 / (K - 1)^2)
  O <- table(factor(a, 0:(K - 1)), factor(b, 0:(K - 1))) / length(a)
  E <- outer(rowSums(O), colSums(O))
  1 - sum(W * O) / sum(W * E)
}
cat_of <- function(x) cut(x, c(-1, 0, 5, 24), labels = c("0分", "1–5分", "6分以上"))   # 常见的临床分组（示例）
c(`Spearman 相关` = cor(manual, model, method = "spearman"), `二次加权kappa` = wkappa(manual, model),
  `分组一致率` = mean(cat_of(manual) == cat_of(model)))
table(人工 = cat_of(manual), 模型 = cat_of(model))
```

```text
[1] 9
Spearman 相关 二次加权kappa    分组一致率 
    0.8072636     0.8029155     0.7200000 
         模型
人工    0分 1–5分 6分以上
  0分     129    78       0
  1–5分    24   209      31
  6分以上   0     7      22
```

读法：这个模拟模型的 Spearman 相关 0.81，二次加权 kappa 0.80，三个分组的一致率 72%。混淆表显示"0 分"里有 78 人被模型判成"1 到 5 分"。**排序和加权 kappa 都不错的模型，在分组边界处仍然会有不少误分**，所以用于分组分析前，要在边界附近单独评估，并且报告人工评分者之间的一致性作为上限参考。

### 小结

- 用 GWAS Catalog 的 API 核对编号和出处，下载时先看文件头、尺度和基因组版本。
- 检索量可以量化"没人做过"，但命中条数必须人工读；只能说没找到。
- Open Targets 的可成药性只说明证据成熟程度，不说明疗效。
- AAC 评分模型要同时看相关、加权 kappa 和分组一致率。

---

[← 上一篇：09 · 个体级队列分析：CHIP × 骨密度 × 骨折，以及共病](09-ukb-cohort-and-comorbidity.md) · [目录](README.md) · [下一篇：11 · 把它们串起来：三阶段流水线与“三取二”立项规则 →](11-pipeline-and-go-nogo.md)
