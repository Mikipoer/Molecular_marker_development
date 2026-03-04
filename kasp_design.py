      
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
KASP 设计一体化（跨平台 + Windows 安全临时文件名 版）
- 双等位过滤（ALT 不含逗号；REF/ALT ∈ {A,C,G,T}，长度=1）
- 预检：染色体存在、坐标有效、参考碱基与 REF 一致、参考不为 N、ALT≠REF
- 设计：固定 20bp 左引物（SNP 为 3' 末端），在下游窗口枚举右引物，筛 Tm/GC/产物（≤60bp）
- 特异性检查（可选）：提供 --db 时启用
  * F1/F2 批量 BLAST（一次）
  * R 候选批量 BLAST（一次）
  * F↔R 命中配对做 in-silico PCR，统计 on/off-target
- 多方案：每个位点可输出前 K 套（--topk）
  排序：off-target 总数 → R 高相似命中数 → 理化评分(|Tm-60|+|GC-0.5|)
- 输出
  * kasp_primers.csv：主表（含加/不加接头的 F1/F2、R、Tm/GC、产物、命中统计、Design_Rank）
  * kasp_blast_hits.tsv：命中明细（R 行含 design_rank）
  * kasp_amplicons.tsv：扩增子明细（含 design_rank）
  * snp_check.tsv / snp_pass.tsv / snp_filtered_out.tsv（如有）
  * kasp_locus_log.tsv：逐位点日志（候选数、失败原因、最佳指标等，周期性落盘）
  * kasp_run.log：运行日志（含进度）

可作为脚本：
  python kasp_design.py -f Zm-B73-REFERENCE-GRAMENE-4.0.fasta -s snp.tsv \
    --db B73db --threads 8 -k 3 \
    -o kasp_primers.csv --hits kasp_blast_hits.tsv --amplicons kasp_amplicons.tsv \
    --locus-log kasp_locus_log.tsv --run-log kasp_run.log --progress-every 10

或作为库（Python 内调用）：
  from kasp_design import run_kasp
  run_kasp(fasta, snp, dbprefix="B73db", threads=8, topk=3)
"""

from __future__ import annotations
import os, sys, re, shutil, subprocess, tempfile, logging
from typing import Dict, Optional, Tuple, List
import pandas as pd
from pyfaidx import Fasta

# ========== primer3 Tm：优先用新接口 calc_tm，旧版回退 ==========
try:
    from primer3.bindings import calc_tm as _calc_tm
except Exception:
    from primer3.bindings import calcTm as _calc_tm

# ================== 可调参数（核心） ==================
FAM_TAIL = "GAAGGTGACCAAGTTCATGCT"   # F1 5' 尾
HEX_TAIL = "GAAGGTCGGAGTCAACGGATT"  # F2 5' 尾

LEFT_LEN = 20                        # 左引物固定长度（SNP 为第 20 位）
RIGHT_LEN_RANGE = (18, 25)           # 右引物长度范围
TM_RANGE = (58.0, 62.0)              # Tm（℃）
GC_RANGE = (0.40, 0.60)              # GC 比例
PRODUCT_MAX = 60                     # 产物长度上限（bp）
RIGHT_SCAN_MAX = 40                  # 下游扫描窗口（bp）

MM_GAP_THRESHOLD = 5                 # 高相似阈值（mismatch+gaps ≤ 5）
DEFAULT_THREADS = 8
PROGRESS_EVERY = 20                  # 每处理 N 个位点打印进度并落盘位点日志

# BLAST 小优化（可按需调整）
BLAST_PERC_ID = 75                   # 最低百分比身份
BLAST_MAX_TARGET_SEQS = 2000         # 防止重复区输出爆炸
BLAST_MAX_HSPS = 1                   # 每命中只保留一个 HSP

CHR_PREFIX_TRIALS = ["", "chr", "Chr", "chr0", "Chr0"]
VALID_BASE = set("ACGT")

# ================== 日志 ==================
def setup_logger(log_path: Optional[str]) -> logging.Logger:
    logger = logging.getLogger("KASP")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")

    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    logger.addHandler(sh)

    if log_path:
        fh = logging.FileHandler(log_path, mode="w", encoding="utf-8")
        fh.setFormatter(fmt)
        logger.addHandler(fh)
    return logger

# ================== 工具函数 ==================
def revcomp(seq: str) -> str:
    t = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(t)[::-1]

def gc(seq: str) -> float:
    s = seq.upper()
    return (s.count("G") + s.count("C")) / max(1, len(s))

def tm_value(seq: str) -> float:
    return _calc_tm(seq)

def find_chrom(fa: Fasta, chrom_raw: str) -> str:
    if chrom_raw in fa: return chrom_raw
    for p in CHR_PREFIX_TRIALS:
        cand = f"{p}{chrom_raw}"
        if cand in fa: return cand
    raise KeyError(f"FASTA 中找不到染色体：{chrom_raw}")

def fa_base_1b(fa: Fasta, chrom: str, pos1: int) -> str:
    return str(fa[chrom][pos1-1:pos1]).upper()

def fa_window_1b(fa: Fasta, chrom: str, start1: int, end1: int) -> str:
    start1 = max(1, start1)
    return str(fa[chrom][start1-1:end1]).upper()

# ================== 读取 + 双等位过滤 ==================
def read_biallelic_snp_tsv(snp_tsv: str) -> pd.DataFrame:
    """
    仅按空白分隔读取 snp.tsv（表头需含 CHR POS REF ALT）；
    只保留双等位：ALT 不含逗号，且 REF/ALT ∈ {A,C,G,T} 且长度=1。
    """
    df = pd.read_csv(snp_tsv, sep=r"\s+", engine="python", comment="#", header=0)
    df = df.iloc[:, :4].copy()
    df.columns = ["CHR", "POS", "REF", "ALT"]

    df["CHR"] = df["CHR"].astype(str)
    df["POS"] = df["POS"].astype(int)
    df["REF"] = df["REF"].astype(str).str.upper().str.strip()
    df["ALT"] = df["ALT"].astype(str).str.upper().str.strip()

    mask = (
        (~df["ALT"].str.contains(",", regex=False))
        & (df["REF"].str.len() == 1)
        & (df["ALT"].str.len() == 1)
        & (df["REF"].isin(VALID_BASE))
        & (df["ALT"].isin(VALID_BASE))
    )
    kept = df[mask].copy()
    bad = df.loc[~mask].copy()
    if len(bad):
        bad.to_csv("snp_filtered_out.tsv", sep="\t", index=False)
    return kept

# ================== 预检 ==================
def precheck_all(
    fasta: str,
    snp_tsv: str,
    out_report: str,
    out_pass: str,
    logger: logging.Logger,
    locus_log: dict,
) -> pd.DataFrame:
    """返回仅 PASS 的位点（CHR POS REF ALT），并写 snp_check.tsv / snp_pass.tsv。"""
    fa = Fasta(fasta, as_raw=True, sequence_always_upper=True)
    snp = read_biallelic_snp_tsv(snp_tsv)

    rows = []
    for _, r in snp.iterrows():
        c_raw, pos, REF, ALT = r["CHR"], int(r["POS"]), r["REF"], r["ALT"]
        status, reasons = "PASS", []
        chrom_in, ref_on_ref = "", ""

        try:
            chrom_in = find_chrom(fa, c_raw)
        except KeyError:
            status, chrom_in = "FAIL", ""
            reasons.append("chrom_not_found")

        if status == "PASS":
            try:
                ref_on_ref = fa_base_1b(fa, chrom_in, pos)
            except Exception as e:
                status = "FAIL"
                reasons.append(f"pos_out_of_range:{e}")

        if status == "PASS":
            if ref_on_ref == "N":
                status = "FAIL"; reasons.append("reference_is_N")
            if REF == ALT:
                status = "FAIL"; reasons.append("ALT_equals_REF")
            if ref_on_ref != REF:
                status = "FAIL"; reasons.append(f"REF_mismatch(ref={ref_on_ref})")

        rows.append({
            "CHR": c_raw, "POS": pos,
            "REF_input": REF, "ALT_input": ALT,
            "Chrom_in_FASTA": chrom_in, "Ref_on_FASTA": ref_on_ref,
            "Status": status, "Reason": ";".join(reasons)
        })

        key = (str(c_raw), int(pos))
        locus_log[key] = {
            "CHR": c_raw, "POS": pos,
            "precheck_status": status,
            "precheck_reason": ";".join(reasons) if reasons else "",
            "left_start": "", "up_len": "", "down_len": "",
            "F1": "", "F2": "", "F1_Tm": "", "F2_Tm": "", "F1_GC": "", "F2_GC": "",
            "cands_total": 0, "cands_physchem": 0,
            "selected_count": 0,
            "best_Product_Size": "", "best_R_hits_mm<=5": "", "best_OffTarget_Total": "",
            "design_status": "SKIPPED", "design_reason": "",
        }

    rep = pd.DataFrame(rows).sort_values(["CHR","POS"])
    rep.to_csv(out_report, sep="\t", index=False)

        # 为避免 Pandas FutureWarning，改用 merge：
    passed = rep[rep["Status"]=="PASS"][["CHR","POS"]].merge(
        snp[["CHR","POS","REF","ALT"]], on=["CHR","POS"], how="left"
    )
    passed.to_csv(out_pass, sep="\t", index=False)

    logger.info(f"[预检] PASS {len(passed)} / {len(rep)}，报告：{out_report}，可用位点：{out_pass}")
    return passed

# ================== BLAST（批量） ==================
_COLS = ["qseqid","sseqid","sstart","send","sstrand","length","pident","nident","mismatch","gaps","evalue","bitscore"]

def _post_blast_df(out_path: str) -> pd.DataFrame:
    if not os.path.exists(out_path) or os.path.getsize(out_path)==0:
        return pd.DataFrame(columns=_COLS)
    df = pd.read_csv(out_path, sep="\t", header=None, names=_COLS)
    df["mm_gap"]    = df["mismatch"].astype(int) + df["gaps"].astype(int)
    df["hit_start"] = df[["sstart","send"]].min(axis=1).astype(int)
    df["hit_end"]   = df[["sstart","send"]].max(axis=1).astype(int)
    return df

# ✅ Windows 安全的临时文件名（仅用于文件名；qseqid 仍保留原样）
def _safe_label_for_filename(label: str) -> str:
    # Windows 禁止 <>:"/\|?* 等字符；统一替换为下划线
    return re.sub(r'[^A-Za-z0-9._-]+', '_', label)

def run_blast_batch(seq_map: Dict[str,str], dbprefix: str, tmpdir: str, qlabel: str, threads: int) -> Dict[str,pd.DataFrame]:
    """
    批量 BLAST：seq_map = {qid: seq, ...}
    返回：{qid: hits_df}（若无命中返回空 DF）
    """
    qfa = os.path.join(tmpdir, f"{_safe_label_for_filename(qlabel)}.fa")
    with open(qfa, "w") as f:
        for qid, s in seq_map.items():
            f.write(f">{qid}\n{s}\n")  # qseqid 保留原始含冒号/竖线，便于后续分组
    out = os.path.join(tmpdir, f"{_safe_label_for_filename(qlabel)}.tsv")
    cmd = [
        "blastn","-task","blastn-short",
        "-query", qfa, "-db", dbprefix,
        "-evalue","1000",
        "-word_size","7",
        "-perc_identity", str(BLAST_PERC_ID),
        "-max_target_seqs", str(BLAST_MAX_TARGET_SEQS),
        "-max_hsps", str(BLAST_MAX_HSPS),
        "-dust","no","-soft_masking","false",
        "-outfmt","6 qseqid sseqid sstart send sstrand length pident nident mismatch gaps evalue bitscore",
        "-num_threads", str(threads),
    ]
    with open(out, "w") as fo:
        subprocess.run(cmd, check=True, stdout=fo)
    all_df = _post_blast_df(out)
    if all_df.empty:
        return {qid: pd.DataFrame(columns=_COLS) for qid in seq_map}
    return { qid: g.copy() for qid, g in all_df.groupby("qseqid") }

# ================== in-silico PCR 配对 ==================
def pair_amplicons(
    f_hits: pd.DataFrame, r_hits: pd.DataFrame,
    target_chr: str, left_start: int, target_pos: int, max_prod: int,
    thr: int = MM_GAP_THRESHOLD
) -> Tuple[List[dict], int, int, Optional[int]]:
    """
    将 F 命中与 R 命中配对形成扩增子：
      - 同染色体；方向 F(+)/R(-) 或 F(-)/R(+)
      - size = r5 - left_start + 1 ≤ max_prod
      - on_target：在目标染色体，覆盖 target_pos
    返回：(amps列表, on_count, off_count, best_on_size)
    """
    if f_hits is None or r_hits is None or f_hits.empty or r_hits.empty:
        return [], 0, 0, None
    f_ok = f_hits[f_hits["mm_gap"]<=thr].copy()
    r_ok = r_hits[r_hits["mm_gap"]<=thr].copy()
    if f_ok.empty or r_ok.empty:
        return [], 0, 0, None

    amps, on, off, best = [], 0, 0, None
    for _, fh in f_ok.iterrows():
        for _, rh in r_ok.iterrows():
            if fh["sseqid"] != rh["sseqid"]:
                continue
            if fh["sstrand"]=="plus" and rh["sstrand"]=="minus":
                r5 = int(rh["hit_end"])    # R 的 5' 端
            elif fh["sstrand"]=="minus" and rh["sstrand"]=="plus":
                r5 = int(rh["hit_start"])
            else:
                continue
            size = r5 - int(left_start) + 1
            if size<=0 or size>max_prod:
                continue
            on_target = (fh["sseqid"]==target_chr and int(left_start)<=target_pos<=r5)
            if on_target:
                on += 1
                if best is None or size < best:
                    best = size
            else:
                off += 1
            amps.append({
                "sseqid": fh["sseqid"], "start": int(left_start), "end": r5, "size": int(size),
                "f_mm_gap": int(fh["mm_gap"]), "r_mm_gap": int(rh["mm_gap"]), "on_target": on_target
            })
    return amps, on, off, best

# ================== 设计（批量 BLAST + 进度日志） ==================
def _is_null(path_like: Optional[str]) -> bool:
    if not path_like:
        return True
    p = str(path_like).strip()
    return p in {"", "-", "NUL", "/dev/null", os.devnull}

def have_blast() -> bool:
    return shutil.which("blastn") is not None

def design_kasp_lt60(
    fasta: str,
    passed_df: pd.DataFrame,
    out_csv: str,
    dbprefix: Optional[str],
    hits_tsv: Optional[str],
    amp_tsv: Optional[str],
    threads: int,
    top_k: int,
    logger: logging.Logger,
    locus_log: dict,
    progress_every: int,
    locus_log_path: Optional[str],
):
    fa = Fasta(fasta, as_raw=True, sequence_always_upper=True)
    use_blast = (dbprefix is not None) and have_blast()

    out_rows, hits_rows, amp_rows = [], [], []
    total, done = len(passed_df), 0

    empty_hits = pd.DataFrame(columns=_COLS)

    with tempfile.TemporaryDirectory() as tmpdir:
        for _, r in passed_df.iterrows():
            c_raw, pos, REF, ALT = r["CHR"], int(r["POS"]), r["REF"].upper(), r["ALT"].upper()
            key = (str(c_raw), int(pos))

            # 染色体 + REF 二次核对
            try:
                chrom = find_chrom(fa, c_raw)
            except KeyError:
                if key in locus_log: locus_log[key]["design_reason"] = "chrom_not_found_in_design"
                done += 1; _maybe_progress(logger, done, total, progress_every, locus_log, locus_log_path)
                continue
            if fa_base_1b(fa, chrom, pos) != REF:
                if key in locus_log: locus_log[key]["design_reason"] = "REF_mismatch_in_design"
                done += 1; _maybe_progress(logger, done, total, progress_every, locus_log, locus_log_path)
                continue

            # 左引物
            left_start = pos - (LEFT_LEN - 1)
            if left_start < 1:
                if key in locus_log: locus_log[key]["design_reason"] = "UPSTREAM_TOO_SHORT"
                done += 1; _maybe_progress(logger, done, total, progress_every, locus_log, locus_log_path)
                continue
            up = fa_window_1b(fa, chrom, left_start, pos-1)
            if len(up) != (LEFT_LEN - 1):
                if key in locus_log: locus_log[key]["design_reason"] = "UPSTREAM_TRUNCATED"
                done += 1; _maybe_progress(logger, done, total, progress_every, locus_log, locus_log_path)
                continue

            f1 = up + REF
            f2 = up + ALT
            f1_tm, f2_tm = tm_value(f1), tm_value(f2)
            f1_gc, f2_gc = gc(f1), gc(f2)

            down = fa_window_1b(fa, chrom, pos+1, pos+RIGHT_SCAN_MAX)

            # 右引物候选（局部索引）
            cands: List[Tuple[str,int,int,float,float,int,float]] = []
            cands_total, cands_physchem = 0, 0
            for L in range(RIGHT_LEN_RANGE[0], RIGHT_LEN_RANGE[1]+1):
                for end_off in range(1, len(down)+1):  # 1-based 局部结束位置
                    start_off = end_off - L + 1
                    if start_off < 1:
                        continue
                    cands_total += 1
                    seg = down[start_off-1:end_off]
                    rseq = revcomp(seg)
                    rtm, rgc = tm_value(rseq), gc(rseq)
                    if not (TM_RANGE[0] <= rtm <= TM_RANGE[1] and GC_RANGE[0] <= rgc <= GC_RANGE[1]):
                        continue
                    r5_abs = pos + end_off
                    prod = r5_abs - left_start + 1       # = end_off + LEFT_LEN
                    if prod <= 0 or prod > PRODUCT_MAX:
                        continue
                    score = abs(60-rtm) + abs(0.5-rgc)
                    cands.append((rseq, r5_abs, L, rtm, rgc, prod, score))
                    cands_physchem += 1

            if key in locus_log:
                locus_log[key].update({
                    "left_start": left_start, "up_len": len(up), "down_len": len(down),
                    "F1": f1, "F2": f2,
                    "F1_Tm": round(f1_tm,2), "F2_Tm": round(f2_tm,2),
                    "F1_GC": round(f1_gc,3), "F2_GC": round(f2_gc,3),
                    "cands_total": cands_total, "cands_physchem": cands_physchem
                })

            if not cands:
                if key in locus_log: locus_log[key]["design_reason"] = "NO_R_CANDIDATE_PHYS"
                done += 1; _maybe_progress(logger, done, total, progress_every, locus_log, locus_log_path)
                continue

            cands.sort(key=lambda x: x[-1])  # 先按理化评分粗排
            locus = f"{c_raw}:{pos}"

            # ===== 特异性（批量 BLAST） =====
            if use_blast:
                # F 批量
                f_map = {f"{locus}|F1": f1, f"{locus}|F2": f2}
                f_hits_map = run_blast_batch(f_map, dbprefix, tmpdir, f"{locus}|Fbatch", threads)
                f1_hits = f_hits_map.get(f"{locus}|F1", empty_hits)
                f2_hits = f_hits_map.get(f"{locus}|F2", empty_hits)

                F1_cnt = int((f1_hits["mm_gap"]<=MM_GAP_THRESHOLD).sum()) if not f1_hits.empty else 0
                F2_cnt = int((f2_hits["mm_gap"]<=MM_GAP_THRESHOLD).sum()) if not f2_hits.empty else 0

                if not f1_hits.empty:
                    tmp = f1_hits.copy(); tmp.insert(1,"primer","F1"); tmp.insert(2,"seq",f1); tmp.insert(3,"design_rank","-")
                    hits_rows.append(tmp)
                if not f2_hits.empty:
                    tmp = f2_hits.copy(); tmp.insert(1,"primer","F2"); tmp.insert(2,"seq",f2); tmp.insert(3,"design_rank","-")
                    hits_rows.append(tmp)

                # R 批量
                r_list = [(f"{locus}|R{idx}", rseq, r5, L, rtm, rgc, prod, score)
                          for idx, (rseq, r5, L, rtm, rgc, prod, score) in enumerate(cands, start=1)]
                r_seq_map = {qid: rseq for (qid, rseq, *_rest) in r_list}
                r_hits_map = run_blast_batch(r_seq_map, dbprefix, tmpdir, f"{locus}|Rbatch", threads)

                selected = []
                for (qid, rseq, r5, L, rtm, rgc, prod, score) in r_list:
                    r_hits = r_hits_map.get(qid, empty_hits)
                    amp1,on1,off1,_ = pair_amplicons(f1_hits, r_hits, chrom, left_start, pos, PRODUCT_MAX)
                    amp2,on2,off2,_ = pair_amplicons(f2_hits, r_hits, chrom, left_start, pos, PRODUCT_MAX)
                    R_cnt = int((r_hits["mm_gap"]<=MM_GAP_THRESHOLD).sum()) if not r_hits.empty else 0
                    if on1>=1 and on2>=1:
                        key_metrics = (off1+off2, R_cnt, score)
                        selected.append({
                            "rseq": rseq, "rtm": rtm, "rgc": rgc, "prod": prod,
                            "key": key_metrics, "r_hits": r_hits,
                            "amp1": amp1, "amp2": amp2,
                            "off_total": off1+off2, "R_cnt": R_cnt
                        })
                if not selected:
                    if key in locus_log: locus_log[key]["design_reason"] = "NO_ON_TARGET_WITH_BLAST"
                    done += 1; _maybe_progress(logger, done, total, progress_every, locus_log, locus_log_path)
                    continue

                selected.sort(key=lambda d: d["key"])
                seen, top = set(), []
                for d in selected:
                    if d["rseq"] in seen: continue
                    top.append(d); seen.add(d["rseq"])
                    if len(top) >= top_k: break

                for rank, d in enumerate(top, start=1):
                    out_rows.append({
                        "Design_Rank": rank,
                        "CHR": c_raw, "POS": pos, "REF": REF, "ALT": ALT,
                        "F1_noTail": f1, "F2_noTail": f2, "R": d["rseq"],
                        "F1_withTail": FAM_TAIL + f1, "F2_withTail": HEX_TAIL + f2,
                        "F1_Tm": round(f1_tm,2), "F2_Tm": round(f2_tm,2), "R_Tm": round(d["rtm"],2),
                        "F1_GC": round(f1_gc,3), "F2_GC": round(f2_gc,3), "R_GC": round(d["rgc"],3),
                        "Product_Size": int(d["prod"]),
                        "F1_hits_mm<=5": F1_cnt, "F2_hits_mm<=5": F2_cnt,
                        "R_hits_mm<=5": d["R_cnt"], "OffTarget_Total": d["off_total"]
                    })
                    if d["r_hits"] is not None and not d["r_hits"].empty:
                        rh = d["r_hits"].copy()
                        rh.insert(1,"primer","R"); rh.insert(2,"seq",d["rseq"]); rh.insert(3,"design_rank",rank)
                        hits_rows.append(rh)
                    for a in d["amp1"]:
                        amp_rows.append({"locus": locus, "design_rank": rank, "pair": "F1-R", **a})
                    for a in d["amp2"]:
                        amp_rows.append({"locus": locus, "design_rank": rank, "pair": "F2-R", **a})

                if key in locus_log:
                    locus_log[key].update({
                        "selected_count": len(top),
                        "best_Product_Size": top[0]["prod"] if top else "",
                        "best_R_hits_mm<=5": top[0]["R_cnt"] if top else "",
                        "best_OffTarget_Total": top[0]["off_total"] if top else "",
                        "design_status": "DESIGNED", "design_reason": ""
                    })

            else:
                # 不做 BLAST：仅按理化评分取前 K
                seen, top = set(), []
                for (rseq, r5, L, rtm, rgc, prod, score) in cands:
                    if rseq in seen: continue
                    top.append((rseq, rtm, rgc, prod, score)); seen.add(rseq)
                    if len(top) >= top_k: break
                for rank, (rseq, rtm, rgc, prod, score) in enumerate(top, start=1):
                    out_rows.append({
                        "Design_Rank": rank,
                        "CHR": c_raw, "POS": pos, "REF": REF, "ALT": ALT,
                        "F1_noTail": f1, "F2_noTail": f2, "R": rseq,
                        "F1_withTail": FAM_TAIL + f1, "F2_withTail": HEX_TAIL + f2,
                        "F1_Tm": round(f1_tm,2), "F2_Tm": round(f2_tm,2), "R_Tm": round(rtm,2),
                        "F1_GC": round(f1_gc,3), "F2_GC": round(f2_gc,3), "R_GC": round(rgc,3),
                        "Product_Size": int(prod),
                        "F1_hits_mm<=5": 0, "F2_hits_mm<=5": 0,
                        "R_hits_mm<=5": "", "OffTarget_Total": ""
                    })
                if key in locus_log:
                    locus_log[key].update({
                        "selected_count": len(top),
                        "best_Product_Size": top[0][3] if top else "",
                        "design_status": "DESIGNED" if top else "SKIPPED",
                        "design_reason": "" if top else "NO_R_SELECTED"
                    })

            done += 1
            _maybe_progress(logger, done, total, progress_every, locus_log, locus_log_path)

    # ---- 写主结果 ----
    if not out_rows:
        print("没有成功设计的位点。", file=sys.stderr)
        sys.exit(2)

    out = pd.DataFrame(out_rows)
    cols_order = [
        "Design_Rank","CHR","POS","REF","ALT",
        "F1_noTail","F2_noTail","R",
        "F1_withTail","F2_withTail",
        "F1_Tm","F2_Tm","R_Tm","F1_GC","F2_GC","R_GC",
        "Product_Size","F1_hits_mm<=5","F2_hits_mm<=5","R_hits_mm<=5","OffTarget_Total"
    ]
    out = out[cols_order]
    out.to_csv(out_csv, index=False)

    # ---- 命中明细 ----
    if not _is_null(hits_tsv):
        if hits_rows:
            hits = pd.concat(hits_rows, ignore_index=True)
            hits.to_csv(hits_tsv, sep="\t", index=False)
        else:
            open(hits_tsv, "w").close()

    # ---- 扩增子明细 ----
    if not _is_null(amp_tsv):
        if amp_rows:
            amp = pd.DataFrame(amp_rows)
            amp.to_csv(amp_tsv, sep="\t", index=False)
        else:
            open(amp_tsv, "w").close()

    if locus_log_path:
        # 兜底：若有 SKIPPED 但无原因，补 EARLY_SKIP_UNKNOWN
        for v in locus_log.values():
            if v.get("design_status")=="SKIPPED" and not v.get("design_reason"):
                v["design_reason"] = "EARLY_SKIP_UNKNOWN"
        write_locus_log(locus_log, locus_log_path, logger)

    logger.info(f"[完成] 主结果: {out_csv}；命中: {hits_tsv}；扩增子: {amp_tsv}")


def _maybe_progress(
    logger: logging.Logger, done: int, total: int, progress_every: int,
    locus_log: dict, locus_log_path: Optional[str]
):
    if progress_every <= 0: return
    if (done % progress_every == 0) or (done == total):
        logger.info(f"[进度] {done}/{total} loci processed")
        if locus_log_path:
            write_locus_log(locus_log, locus_log_path, logger=None)

# ================== 逐位点日志 ==================
def write_locus_log(locus_log: dict, path: str, logger: Optional[logging.Logger]):
    if not path: return
    rows = [v for _, v in locus_log.items()]
    if not rows:
        open(path, "w").close(); return
    df = pd.DataFrame(rows)
    cols = [
        "CHR","POS",
        "precheck_status","precheck_reason",
        "left_start","up_len","down_len",
        "F1","F2","F1_Tm","F2_Tm","F1_GC","F2_GC",
        "cands_total","cands_physchem",
        "selected_count","best_Product_Size","best_R_hits_mm<=5","best_OffTarget_Total",
        "design_status","design_reason",
    ]
    cols = [c for c in cols if c in df.columns]
    df = df[cols].sort_values(["CHR","POS"])
    df.to_csv(path, sep="\t", index=False)
    if logger: logger.info(f"[日志] 位点级日志刷新：{path}")

# ================== 可选：自动建库 ==================
def ensure_blastdb_if_needed(fasta: str, dbprefix: str) -> None:
    """
    如果 dbprefix.* 不存在，自动调用 makeblastdb 建库。
    注意：需要系统 PATH 里有 makeblastdb（Windows 也支持）。
    """
    if os.path.exists(dbprefix + ".nsq"):
        return
    cmd = ["makeblastdb", "-in", fasta, "-dbtype", "nucl", "-parse_seqids", "-out", dbprefix]
    print("[info] building BLAST DB:", " ".join(cmd))
    subprocess.run(cmd, check=True)

# ================== Python API ==================
def run_kasp(
    fasta: str,
    snp_tsv: str,
    dbprefix: Optional[str] = None,   # 传 None 则跳过 BLAST（更快）
    threads: int = DEFAULT_THREADS,
    topk: int = 1,
    out_csv: str = "kasp_primers.csv",
    hits_tsv: Optional[str] = "kasp_blast_hits.tsv",
    amp_tsv: Optional[str] = "kasp_amplicons.tsv",
    locus_log_path: str = "kasp_locus_log.tsv",
    run_log_path: str = "kasp_run.log",
    progress_every: int = PROGRESS_EVERY,
    auto_make_db: bool = False,       # True 时，在 Python 内自动 makeblastdb（若缺失）
):
    logger = setup_logger(run_log_path)
    logger.info("====== KASP 设计启动（Python API）======")
    logger.info(f"FASTA: {fasta}")
    logger.info(f"SNP:   {snp_tsv}")
    if dbprefix:
        logger.info(f"BLAST DB: {os.path.basename(dbprefix)}（threads={threads}）")
        if auto_make_db:
            ensure_blastdb_if_needed(fasta, dbprefix)
    else:
        logger.info("BLAST: 跳过（dbprefix=None）")
    logger.info(f"每个位点返回前 {topk} 套；进度每 {progress_every} 条刷新")

    locus_log: dict = {}
    report, pass_out = "snp_check.tsv", "snp_pass.tsv"
    passed = precheck_all(fasta, snp_tsv, report, pass_out, logger, locus_log)

    design_kasp_lt60(
        fasta=fasta,
        passed_df=passed,
        out_csv=out_csv,
        dbprefix=dbprefix,
        hits_tsv=hits_tsv,
        amp_tsv=amp_tsv,
        threads=threads,
        top_k=topk,
        logger=logger,
        locus_log=locus_log,
        progress_every=progress_every,
        locus_log_path=locus_log_path,
    )
    write_locus_log(locus_log, locus_log_path, logger)
    logger.info("====== KASP 设计完成（Python API）======")

# ================== CLI ==================
def main():
    import argparse
    ap = argparse.ArgumentParser(
        description="KASP 一体化：双等位过滤 + 预检 + 设计（≤60bp）+ 批量 BLAST + 多方案 + 进度/位点日志"
    )
    ap.add_argument("-f","--fasta", required=True, help="参考基因组 FASTA")
    ap.add_argument("-s","--snp",   required=True, help="SNP 表（空白分隔，含表头 CHR POS REF ALT，1-based）")

    ap.add_argument("--report",   default="snp_check.tsv", help="预检报告 TSV")
    ap.add_argument("--pass-out", default="snp_pass.tsv",  help="仅 PASS 位点 TSV")
    ap.add_argument("--check-only", action="store_true",   help="只做预检后退出")

    ap.add_argument("-o","--out", default="kasp_primers.csv",     help="主结果 CSV")
    ap.add_argument("--hits",     default="kasp_blast_hits.tsv",  help="BLAST 命中 TSV（或传 NUL/os.devnull 跳过）")
    ap.add_argument("--amplicons",default="kasp_amplicons.tsv",   help="扩增子 TSV（或传 NUL/os.devnull 跳过）")

    ap.add_argument("--db",       default=None, help="BLAST 数据库前缀（如 B73db；留空跳过 BLAST）")
    ap.add_argument("--threads",  type=int, default=DEFAULT_THREADS, help="BLAST 线程数")
    ap.add_argument("-k","--topk",type=int, default=1, help="每个位点返回的方案数（默认 1）")

    ap.add_argument("--locus-log", default="kasp_locus_log.tsv", help="逐位点日志 TSV")
    ap.add_argument("--run-log",   default="kasp_run.log",       help="运行日志（INFO）")
    ap.add_argument("--progress-every", type=int, default=PROGRESS_EVERY, help="每处理 N 个位点打印进度并刷新位点日志")

    ap.add_argument("--mkdb", action="store_true", help="若 DB 缺失则自动 makeblastdb（需 makeblastdb 在 PATH 中）")

    args = ap.parse_args()

    logger = setup_logger(args.run_log if args.run_log else None)
    logger.info("====== KASP 设计启动 ======")
    logger.info(f"FASTA: {args.fasta}")
    logger.info(f"SNP:   {args.snp}")
    if args.db: logger.info(f"BLAST DB: {args.db}（threads={args.threads}）")
    else:       logger.info("BLAST: 跳过（未提供 --db）")
    logger.info(f"每个位点返回前 {args.topk} 套；进度每 {args.progress_every} 条刷新")

    locus_log: dict = {}

    passed = precheck_all(args.fasta, args.snp, args.report, args.pass_out, logger, locus_log)
    if args.check_only:
        write_locus_log(locus_log, args.locus_log, logger)
        logger.info("仅预检完成（--check-only）。")
        return

    if args.db and args.mkdb:
        ensure_blastdb_if_needed(args.fasta, args.db)

    design_kasp_lt60(
        fasta=args.fasta,
        passed_df=passed,
        out_csv=args.out,
        dbprefix=args.db,
        hits_tsv=args.hits,
        amp_tsv=args.amplicons,
        threads=args.threads,
        top_k=args.topk,
        logger=logger,
        locus_log=locus_log,
        progress_every=args.progress_every,
        locus_log_path=args.locus_log
    )
    logger.info("====== KASP 设计完成 ======")

if __name__ == "__main__":
    main()

    
