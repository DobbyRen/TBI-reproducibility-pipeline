# 基因集来源（data/gmt/）

- 数据库：MSigDB v7.5.1（Human, symbols）
- 下载地址：https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.5.1/msigdb.v7.5.1.symbols.gmt
- 本地文件 md5：3215a7edac29f94ecf9dca4c34427d21
- 拆分时间：2026-10-21 23:39:09
- Hallmark 50 个集合 / KEGG 186 个 / GO_BP 7658 个

说明：MSigDB 官网下载需登录；Broad 开放镜像无需登录且为同一份数据。
所有非人数据集在 L6 之前已通过 09_ortholog.R 映射到人类基因空间，
因此只用人类基因集即可覆盖全部数据集。
