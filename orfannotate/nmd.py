# This currently relies on the 50nt rule
def predict_nmd(orf_end, junctions, strand):
    if not junctions:
        return 'FALSE'

    def _start(f):
        return f.start if hasattr(f, 'start') else f[0]

    def _end(f):
        return f.end if hasattr(f, 'end') else f[1]

    if strand == '+':
        last_donor = _end(junctions[-1])
        distance = last_donor - orf_end
    else:
        first_acceptor = _start(junctions[0])
        distance = orf_end - first_acceptor

    return 'TRUE' if distance > 50 else 'FALSE'