def predict_nmd(orf_end, junctions, strand):
    if not junctions:
        return 'FALSE'
    if strand == '+':
        last_junction = junctions[-1][0]
        distance = orf_end - last_junction
    else:
        last_junction = junctions[0][1]
        distance = last_junction - orf_end
    return 'TRUE' if distance < -50 else 'FALSE'