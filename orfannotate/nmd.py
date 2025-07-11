# This currently relies on the 50nt rule
def predict_nmd(orf_end, junctions, strand):
    if not junctions:
        return 'FALSE'

    if strand == '+':
        last_donor = junctions[-1][0]
        distance = last_donor - orf_end
    else:
        first_acceptor = junctions[0][1]
        distance = orf_end - first_acceptor

    return 'TRUE' if distance > 50 else 'FALSE'