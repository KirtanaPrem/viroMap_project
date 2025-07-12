# -------- Advanced Codon Metrics --------
import numpy as np
from Bio.SeqUtils.CodonUsage import CodonAdaptationIndex

def calc_cai(seq):
    cai_calculator = CodonAdaptationIndex()
    # Load a generic highly‑expressed reference table (E. coli). You can replace with host‑specific table.
    cai_calculator.generate_index()
    return round(cai_calculator.cai_for_gene(seq), 4)

def calc_enc(seq):
    codons = [seq[i:i+3] for i in range(0, len(seq)-2, 3)]
    codon_table = {}
    for c in codons:
        aa = SeqIO.Seq(c).translate(table=1)
        codon_table.setdefault(aa, []).append(c)
    F = []
    for aa, cds in codon_table.items():
        counts = [codons.count(c) for c in cds]
        n = sum(counts)
        if n == 0: continue
        denom = n*(n-1)
        F_aa = (sum(c**2 for c in counts) - n) / denom if denom else 0
        F.append(F_aa)
    enc = 2 + 9/np.mean(F[:2]) + 1/np.mean(F[2:]) if F else np.nan
    return round(enc, 2)

def rare_codon_freq(seq, rare_list={"ATA","AGG","AGA","CGG","CTA","TTA"}):
    codons = [seq[i:i+3] for i in range(0, len(seq)-2, 3)]
    rare = sum(1 for c in codons if c in rare_list)
    return round(rare/len(codons), 4)

def codon_pair_bias(seq):
    pairs = [seq[i:i+6] for i in range(0, len(seq)-5, 3)]
    expected = 1/4096  # uniform expectation (64 codons)^2
    obs = pd.Series(pairs).value_counts(normalize=True)
    return round(((obs-expected)**2/expected).sum(), 4)

def volatility(seq):
    import itertools
    syn_sites = 0
    vol = 0
    for i in range(0, len(seq)-2, 3):
        codon = seq[i:i+3]
        aa = SeqIO.Seq(codon).translate(table=1)
        for pos in range(3):
            for nt in "ACGT":
                if nt == codon[pos]: continue
                new_codon = codon[:pos]+nt+codon[pos+1:]
                new_aa = SeqIO.Seq(new_codon).translate(table=1)
                syn_sites += 1
                if new_aa != aa:
                    vol += 1
    return round(vol/syn_sites if syn_sites else 0, 4)

# ---- Calculate metrics ----
cai = calc_cai(sequence)
enc = calc_enc(sequence)
rcf = rare_codon_freq(sequence)
cpb = codon_pair_bias(sequence)
vol = volatility(sequence)

metric_df = pd.DataFrame({
    "Metric": ["CAI", "ENC", "Rare‑Codon Freq", "Codon Pair Bias", "Volatility"],
    "Value": [cai, enc, rcf, cpb, vol]
})
st.dataframe(metric_df, use_container_width=True)
