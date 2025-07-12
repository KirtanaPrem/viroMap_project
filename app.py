with tabs[1]:  # Codon Bias Tab
    st.markdown("## Codon Bias Metrics (Real-time Calculations)")

    if fasta:
        def get_codon_list(seq):
            return [seq[i:i+3] for i in range(0, len(seq)-2, 3)]

        def rare_codon_freq(codons, rare_set={"ATA", "AGG", "AGA", "CGG", "CTA", "TTA"}):
            return round(sum(1 for c in codons if c in rare_set) / len(codons), 4)

        def codon_pair_bias(codons):
            pairs = [codons[i]+codons[i+1] for i in range(len(codons)-1)]
            expected = 1 / (64*64)
            from collections import Counter
            counts = Counter(pairs)
            total = sum(counts.values())
            observed = [c/total for c in counts.values()]
            cpb = sum((o - expected) ** 2 / expected for o in observed)
            return round(cpb, 4)

        def volatility(codons):
            transitions = {"A": "G", "G": "A", "C": "T", "T": "C"}
            amino_changes = 0
            total = 0
            for codon in codons:
                if len(codon) != 3:
                    continue
                for i in range(3):
                    for base in "ACGT":
                        if base != codon[i]:
                            mutant = codon[:i] + base + codon[i+1:]
                            total += 1
                            if mutant != codon:
                                amino_changes += 1
            return round(amino_changes / total, 4) if total else 0

        def simple_cai(seq):
            # Simplified CAI using a reference E. coli codon usage table
            ref = {
                "GCT": 0.18, "GCC": 0.55, "GCA": 0.15, "GCG": 0.12,  # Ala
                "CGT": 0.04, "CGC": 0.11, "CGA": 0.11, "CGG": 0.21, "AGA": 0.21, "AGG": 0.32,  # Arg
                # You can expand this for all codons
            }
            codons = get_codon_list(seq)
            weights = [ref.get(codon, 0.01) for codon in codons]
            from math import exp, log
            return round(exp(sum(log(w) for w in weights) / len(weights)), 4) if weights else 0

        def enc(codons):
            from collections import Counter
            total = len(codons)
            freqs = Counter(codons)
            n = sum(freqs.values())
            homo = sum(c * (c - 1) for c in freqs.values())
            return round(61 / (1 + homo / (n * (n - 1))), 2) if n > 1 else 0

        clean_seq = "".join(fasta.splitlines()[1:]).replace(" ", "").replace("\n", "").upper()
        codons = get_codon_list(clean_seq)

        metrics = {
            "Codon Adaptation Index (CAI)": simple_cai(clean_seq),
            "Rare Codon Frequency": rare_codon_freq(codons),
            "Codon Pair Bias (CPB)": codon_pair_bias(codons),
            "Codon Volatility": volatility(codons),
            "Effective Number of Codons (ENC)": enc(codons)
        }

        st.dataframe(pd.DataFrame(metrics.items(), columns=["Metric", "Value"]))
    else:
        st.warning("Please search and select a virus to view codon bias metrics.")
