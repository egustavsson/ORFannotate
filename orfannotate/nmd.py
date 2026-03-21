# This currently relies on the 50nt rule
def predict_nmd(orf_end, junctions, strand):

    # Need at least 2 exons to have an exon-exon junction
    if len(junctions) < 2:
        return 'FALSE'

    def _start(f):
        return f.start if hasattr(f, 'start') else f[0]

    def _end(f):
        return f.end if hasattr(f, 'end') else f[1]

    # Evaluates whether the last downstream junction in the 3'UTR is located > 50bp from the last CDS exon:
    if strand == '+':
        ej_positions = [_end(junctions[i]) for i in range(len(junctions) - 1)]
        downstream_ejs = [ej for ej in ej_positions if ej > orf_end]
        if not downstream_ejs:
            return 'FALSE'
        last_ej = downstream_ejs[-1]
        return 'TRUE' if last_ej - orf_end > 50 else 'FALSE'
    else:
        ej_positions = [_start(junctions[i]) for i in range(1, len(junctions))]
        downstream_ejs = [ej for ej in ej_positions if ej < orf_end]
        if not downstream_ejs:
            return 'FALSE'
        last_ej = downstream_ejs[0]
        return 'TRUE' if orf_end - last_ej > 50 else 'FALSE'