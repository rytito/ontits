#!/usr/bin/env python3
"""Demultiplexado de barcodes duales estilo Illumina en lecturas ONT.

Estructura esperada de cada amplicón (hebra sentido):

    [P5 + sitio Read1 (58 nt)] [barcode F (12 nt)] [primer F]
        ... inserto ITS ...
    [rc(primer R)] [rc(barcode R) (12 nt)] [rc(sitio Read2 + P7) (58 nt)]

ONT lee ambas hebras, así que la ronda 1 usa --revcomp de cutadapt para
orientar cada lectura a la vez que la demultiplexa por su barcode frontal.
El barcode solo se acepta ADYACENTE a su primer (barcode+primer como un
único adaptador): con solo 4 barcodes Golay compartidos entre ambos
extremos, un barcode aislado no distingue el extremo F del R -- el primer
adyacente sí.

    ronda 1: -g <barcodeF+primerF>  (orienta, asigna F, recorta el 5')
    ronda 2: -a rc(<barcodeR+primerR>)  (asigna R, recorta el 3')

Una lectura sin ambos barcodes va a unassigned/ -- el equivalente de
--barcode-both-ends de dorado: suprime el index hopping, que con
indexado combinatorio (mismos 4 Golay en ambos extremos) produciría
asignaciones cruzadas silenciosas.

Uso:
    demux_dual_barcodes.py --fastq all.fastq --barcodes barcodes.tsv \
        --fwd-primer TAGAGGAAGTAAAAGTCGTAA --rev-primer TTYRCTRCGTTCTTCATC \
        --outdir demux [--error-rate 0.15] [--min-overlap 25] [--threads 4]
"""

import argparse
import csv
import subprocess
import sys
from pathlib import Path

COMP = str.maketrans("ACGTURYSWKMBDHVNacgturyswkmbdhvn",
                     "TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn")


def revcomp(seq: str) -> str:
    return seq.translate(COMP)[::-1]


def count_reads(path: Path) -> int:
    if not path.exists():
        return 0
    n = 0
    with open(path) as fh:
        for _ in fh:
            n += 1
    return n // 4


def run_cutadapt(args_list, log_path: Path):
    with open(log_path, "w") as log:
        res = subprocess.run(["cutadapt"] + args_list,
                             stdout=log, stderr=subprocess.STDOUT)
    if res.returncode != 0:
        sys.exit(f"ERROR: cutadapt falló (ver {log_path})")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fastq", required=True)
    ap.add_argument("--barcodes", required=True,
                    help="TSV: SampleID, fwd_id, rev_id, fwd_barcode, rev_barcode")
    ap.add_argument("--fwd-primer", required=True)
    ap.add_argument("--rev-primer", required=True)
    ap.add_argument("--outdir", default="demux")
    ap.add_argument("--error-rate", type=float, default=0.15,
                    help="tasa de error de cutadapt; 0.15 tolera ~4 errores "
                         "en un adaptador de 30-33 nt, adecuado para ONT sup")
    ap.add_argument("--min-overlap", type=int, default=25,
                    help="solapamiento mínimo; cerca de la longitud completa "
                         "del adaptador barcode+primer para evitar "
                         "coincidencias parciales espurias")
    ap.add_argument("--min-length", type=int, default=20,
                    help="descarta lecturas más cortas tras el recorte")
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    unassigned = Path("unassigned")
    logdir = Path("logs")
    r1dir = Path("round1")
    r2dir = Path("round2")
    for d in (outdir, unassigned, logdir, r1dir, r2dir):
        d.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------ hoja
    combos = {}      # (fwd_id, rev_id) -> SampleID
    fwd_bc = {}      # fwd_id -> barcode
    rev_bc = {}      # rev_id -> barcode
    with open(args.barcodes) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            combos[(row["fwd_id"], row["rev_id"])] = row["SampleID"]
            fwd_bc.setdefault(row["fwd_id"], row["fwd_barcode"])
            rev_bc.setdefault(row["rev_id"], row["rev_barcode"])
            if fwd_bc[row["fwd_id"]] != row["fwd_barcode"] \
               or rev_bc[row["rev_id"]] != row["rev_barcode"]:
                sys.exit(f"ERROR: barcode inconsistente para {row['SampleID']}")

    # ------------------------------------------------------ ronda 1: F
    cmd = ["-j", str(args.threads), "-e", str(args.error_rate),
           "-O", str(args.min_overlap), "--revcomp"]
    for fid, bc in sorted(fwd_bc.items()):
        cmd += ["-g", f"{fid}={bc}{args.fwd_primer}"]
    cmd += ["--untrimmed-output", str(unassigned / "no_fwd_barcode.fastq"),
            "-o", str(r1dir / "{name}.fastq"), args.fastq]
    run_cutadapt(cmd, logdir / "round1_fwd.log")

    # ------------------------------------------------------ ronda 2: R
    # La lectura ya está orientada y recortada por 5'; el extremo 3' es
    # inserto + rc(primerR) + rc(barcodeR) + cola rc(Read2+P7). El
    # adaptador rc(barcodeR+primerR) recorta desde el primer hasta el
    # final, cola incluida.
    stats = []
    for fid in sorted(fwd_bc):
        fq = r1dir / f"{fid}.fastq"
        if not fq.exists() or fq.stat().st_size == 0:
            continue
        cmd = ["-j", str(args.threads), "-e", str(args.error_rate),
               "-O", str(args.min_overlap), "-m", str(args.min_length)]
        for rid, bc in sorted(rev_bc.items()):
            cmd += ["-a", f"{rid}={revcomp(bc + args.rev_primer)}"]
        cmd += ["--untrimmed-output",
                str(unassigned / f"no_rev_barcode_{fid}.fastq"),
                "-o", str(r2dir / (fid + "__{name}.fastq")), str(fq)]
        run_cutadapt(cmd, logdir / f"round2_{fid}.log")

        for rid in sorted(rev_bc):
            src = r2dir / f"{fid}__{rid}.fastq"
            n = count_reads(src)
            sample = combos.get((fid, rid))
            if sample:
                # .raw.fastq: el nombre debe diferir del F1R1.fastq que
                # escribe CHOPPER_FILTER -- con stageInMode symlink, un
                # redirect sobre el mismo nombre truncaría este archivo.
                if src.exists():
                    src.rename(outdir / f"{sample}.raw.fastq")
                stats.append((sample, f"{fid}-{rid}", n))
            elif n:
                # combinación no usada en la hoja: index hopping u error
                src.rename(unassigned / f"invalid_{fid}_{rid}.fastq")
                stats.append((f"invalid_{fid}-{rid}", f"{fid}-{rid}", n))

    stats.append(("unassigned_no_fwd", "-",
                  count_reads(unassigned / "no_fwd_barcode.fastq")))
    for fid in sorted(fwd_bc):
        stats.append((f"unassigned_no_rev_{fid}", "-",
                      count_reads(unassigned / f"no_rev_barcode_{fid}.fastq")))

    with open("demux_stats.tsv", "w") as out:
        out.write("sample\tcombo\treads\n")
        for name, combo, n in stats:
            out.write(f"{name}\t{combo}\t{n}\n")

    assigned = sum(n for name, _, n in stats if not name.startswith(("unassigned", "invalid")))
    total = assigned + sum(n for name, _, n in stats if name.startswith(("unassigned", "invalid")))
    pct = 100.0 * assigned / total if total else 0.0
    print(f"demux: {assigned}/{total} lecturas asignadas ({pct:.1f}%)")


if __name__ == "__main__":
    main()
