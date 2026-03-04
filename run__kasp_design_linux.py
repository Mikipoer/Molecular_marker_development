#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import sys

from kasp_design import run_kasp


def main():
    ap = argparse.ArgumentParser(
        description="Run KASP primer design (wrapper of kasp_design.run_kasp) with CLI arguments."
    )

    # 位置参数：你想要的 “python run_design_kasp.py 输入文件 输出文件”
    ap.add_argument("snp_tsv", help="SNP table TSV (whitespace-separated; header: CHR POS REF ALT; 1-based).")
    ap.add_argument("out_csv", help="Output CSV path for primer table (e.g., kasp_primers.csv).")

    # 常用可选参数
    ap.add_argument("-f", "--fasta", required=True, help="Reference genome FASTA (e.g., Zm-B73-REFERENCE...fa).")
    ap.add_argument("--db", default=None, help="BLAST DB prefix (without extensions). If omitted, BLAST is skipped.")
    ap.add_argument("--threads", type=int, default=8, help="Threads for BLAST (default: 8).")
    ap.add_argument("-k", "--topk", type=int, default=1, help="Top K designs per locus (default: 1).")
    ap.add_argument("--progress-every", type=int, default=50, help="Progress log interval (default: 50 loci).")

    # 其它输出文件（默认放到 out_csv 同目录）
    ap.add_argument("--hits-tsv", default=None, help="Optional BLAST hits TSV path. Default: <outdir>/kasp_blast_hits.tsv")
    ap.add_argument("--amp-tsv", default=None, help="Optional amplicons TSV path. Default: <outdir>/kasp_amplicons.tsv")
    ap.add_argument("--locus-log", default=None, help="Optional locus log TSV path. Default: <outdir>/kasp_locus_log.tsv")
    ap.add_argument("--run-log", default=None, help="Optional run log path. Default: <outdir>/kasp_run.log")

    args = ap.parse_args()

    # 统一成绝对路径（避免 cwd 影响）
    fasta = os.path.abspath(args.fasta)
    snp_tsv = os.path.abspath(args.snp_tsv)
    out_csv = os.path.abspath(args.out_csv)

    out_dir = os.path.dirname(out_csv) or os.getcwd()
    os.makedirs(out_dir, exist_ok=True)

    # 给其它输出设默认值（放在 out_csv 同目录）
    hits_tsv = os.path.abspath(args.hits_tsv) if args.hits_tsv else os.path.join(out_dir, "kasp_blast_hits.tsv")
    amp_tsv = os.path.abspath(args.amp_tsv) if args.amp_tsv else os.path.join(out_dir, "kasp_amplicons.tsv")
    locus_log = os.path.abspath(args.locus_log) if args.locus_log else os.path.join(out_dir, "kasp_locus_log.tsv")
    run_log = os.path.abspath(args.run_log) if args.run_log else os.path.join(out_dir, "kasp_run.log")

    dbprefix = os.path.abspath(args.db) if args.db else None

    print("PYTHON:", sys.executable)
    print("FASTA:", fasta)
    print("SNP_TSV:", snp_tsv)
    print("OUT_CSV:", out_csv)
    print("DBPREFIX:", dbprefix if dbprefix else "(BLAST skipped)")
    print("OUT_DIR:", out_dir)

    run_kasp(
        fasta=fasta,
        snp_tsv=snp_tsv,
        dbprefix=dbprefix,  # None 就跳过 BLAST；给前缀就启用
        threads=args.threads,
        topk=args.topk,
        out_csv=out_csv,
        hits_tsv=hits_tsv,
        amp_tsv=amp_tsv,
        locus_log_path=locus_log,
        run_log_path=run_log,
        progress_every=args.progress_every,
        auto_make_db=False,  # 如果你想自动 makeblastdb，可加参数再打开
    )

    print("DONE:", out_csv)


if __name__ == "__main__":
    main()

