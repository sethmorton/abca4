#!/usr/bin/env python3
"""Minimal Marimo app: show a FASTA file in a table."""

import marimo

__generated_with = "0.17.8"
app = marimo.App()


@app.cell
def _():
    import re

    def read_fasta_sequence(path):
        """Read the first sequence from a FASTA file and return it as a single string."""
        seq_lines = []
        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    # header line – skip for now
                    continue
                seq_lines.append(line)
        return "".join(seq_lines)

    def find_motifs(seq):
        """
        Find:
        - ATP-binding P-loop/Walker A–like motifs: GxxxxGK[ST]
        - N-glycosylation motifs: N{P}[ST]{P}
        Returns dict with lists of hits.
        """

        # Walker A / P-loop motif (ATP-binding):
        # GxxxxGKST is canonical; we'll allow S or T in the last position: GxxxxGK[ST]
        atp_pattern = re.compile(r"G....GK[ST]")

        # N-glycosylation motif:
        # N{P}[ST]{P}  -> N-X-S/T, where X != P and last != P
        nglyc_pattern = re.compile(r"N[^P][ST][^P]")

        atp_sites = []
        for m in atp_pattern.finditer(seq):
            start = m.start() + 1  # 1-based indexing
            end = m.end()
            atp_sites.append({
                "motif": m.group(),
                "start": start,
                "end": end,
            })

        nglyc_sites = []
        for m in nglyc_pattern.finditer(seq):
            start = m.start() + 1
            end = m.end()
            nglyc_sites.append({
                "motif": m.group(),
                "start": start,
                "end": end,
            })

        return {"atp_sites": atp_sites, "nglyc_sites": nglyc_sites}

    return read_fasta_sequence, find_motifs


@app.cell
def __(read_fasta_sequence, find_motifs):
    import marimo as mo
    from pathlib import Path

    # Use Path(__file__) to compute the path relative to this notebook
    notebook_dir = Path(__file__).resolve().parent
    fasta_path = notebook_dir.parent / "data_raw" / "sequences" / "ABCA4_P78363.fasta"

    try:
        seq = read_fasta_sequence(str(fasta_path))
        motifs = find_motifs(seq)
        
        atp_output = mo.md("**ATP-binding-like motifs (GxxxxGK[ST]):**\n")
        for site in motifs["atp_sites"]:
            atp_output = atp_output.concat(mo.md(f"- {site['motif']} at {site['start']}-{site['end']}"))
        
        nglyc_output = mo.md("\n**N-glycosylation motifs (N[^P][ST][^P]):**\n")
        for site in motifs["nglyc_sites"]:
            nglyc_output = nglyc_output.concat(mo.md(f"- {site['motif']} at {site['start']}-{site['end']}"))
        
        view = mo.vstack([atp_output, nglyc_output])
    except FileNotFoundError:
        view = mo.md(f"FASTA file not found at `{fasta_path}`")

    return view,

@app.cell
def __(view):
    view
    return

if __name__ == "__main__":
    app.run()
