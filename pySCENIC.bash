#!/bin/bash

#SBATCH --job-name=SCENIC        # 作业名
#SBATCH --partition=64c512g_rocky        # cpu 队列
#SBATCH -n 64                 # 总核数 40
#SBATCH --ntasks-per-node=64   # 每节点核数
#SBATCH --output=%j.out
#SBATCH --error=%j.err

cd /dssg/home/acct-yankepeng/songmei/Project/5.PUMCH_Blood

database1=/dssg/home/acct-yankepeng/songmei/Ref/SCENIC_hg38/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
database2=/dssg/home/acct-yankepeng/songmei/Ref/SCENIC_hg38/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
anno=/dssg/home/acct-yankepeng/songmei/Ref/SCENIC_hg38/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
list=/dssg/home/acct-yankepeng/songmei/Ref/SCENIC_hg38/allTFs_hg38.txt

conda activate pyscenic

pyscenic grn --sparse --method grnboost2 --output sce.adj.csv sample.loom ${list}

pyscenic ctx --output sce.regulons.csv --expression_mtx_fname sample.loom --all_modules --mask_dropouts --mode "dask_multiprocessing" --min_genes 10 --annotations_fname ${anno} sce.adj.csv ${database1} ${database2}

pyscenic aucell --output sce_SCENIC.loom sample.loom sce.regulons.csv

conda deactivate
