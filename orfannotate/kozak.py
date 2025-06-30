# Classify the Kozak sequence around the start codon.
def score_kozak(seq, orf_start):
    if orf_start < 7 or orf_start + 3 >= len(seq):
        return "weak", "NA"

    context = seq[orf_start - 7 : orf_start + 3].upper()  # 6 upstream + start codon + 1
    if len(context) < 10:
        return "weak", context

    base_at_minus3 = context[3]
    base_at_plus4 = context[9]

    if base_at_minus3 in ("A", "G") and base_at_plus4 == "G":
        return "strong", context
    elif base_at_minus3 in ("A", "G") or base_at_plus4 == "G":
        return "moderate", context
    else:
        return "weak", context