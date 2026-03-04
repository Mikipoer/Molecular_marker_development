      
import sys
print("PYTHON:", sys.executable)
from kasp_design import run_kasp
import os

BASE_DIR = r"D:\Desktop\VScode\Molecular_marker_development"  # ★ 修改为你的实际路径
os.chdir(BASE_DIR)

FASTA = os.path.join(BASE_DIR, "Zm-B73-REFERENCE-GRAMENE-4.0.fa")
SNP_TSV = os.path.join(BASE_DIR, "snp_all_fixed_chr.tsv")  # 你实际用的那个SNP文件名

OUT_DIR = os.path.join(BASE_DIR, "kasp_out")
os.makedirs(OUT_DIR, exist_ok=True)

DBPREFIX = os.path.join(BASE_DIR, "B73v4_db")  # ★ BLAST库前缀（不带扩展名）

run_kasp(
    fasta=FASTA,
    snp_tsv=SNP_TSV,
    dbprefix=DBPREFIX,              # ★ 由 None 改为 DBPREFIX
    threads=8,
    topk=5,
    out_csv=os.path.join(OUT_DIR, "kasp_primers_with_blast.csv"),
    hits_tsv=os.path.join(OUT_DIR, "kasp_blast_hits.tsv"),
    amp_tsv=os.path.join(OUT_DIR, "kasp_amplicons.tsv"),
    locus_log_path=os.path.join(OUT_DIR, "kasp_locus_log.tsv"),
    run_log_path=os.path.join(OUT_DIR, "kasp_run.log"),
    progress_every=50
)

print("DONE:", os.path.join(OUT_DIR, "kasp_primers_with_blast.csv"))

    