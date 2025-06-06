from Bio import SeqIO

def generate_summary(transcripts_fasta, orfs, gtf, output_tsv):
    gene_map = {}
    with open(gtf) as fh:
        for line in fh:
            if "\ttranscript\t" in line:
                fields = line.strip().split("\t")
                chrom = fields[0]
                strand = fields[6]
                info = fields[8]
                tid = info.split("transcript_id")[1].split('"')[1]
                gid = info.split("gene_id")[1].split('"')[1]
                gene_map[tid] = (chrom, strand, gid)

    with open(output_tsv, "w") as out:
        header = [
            "transcript_ID", "chrom", "strand", "gene", "coding",
            "ORF_length", "CDS_length", "predicted_NMD",
            "transcript_seq", "CDS_seq", "AA_seq"
        ]
        out.write("\t".join(header) + "\n")
        for record in SeqIO.parse(transcripts_fasta, "fasta"):
            tid = record.id
            tseq = str(record.seq)
            chrom, strand, gid = gene_map.get(tid, ("NA", "NA", "NA"))
            orf = orfs.get(tid)
            if orf:
                cds_seq = orf["nt_seq"]
                aa_seq = orf["aa_seq"]
                out.write("\t".join([
                    tid, chrom, strand, gid, "coding",
                    str(orf["length"]), str(len(cds_seq)), orf["predicted_NMD"],
                    tseq, cds_seq, aa_seq
                ]) + "\n")
            else:
                out.write("\t".join([
                    tid, chrom, strand, gid, "non-coding",
                    "0", "0", "NA", tseq, "", ""
                ]) + "\n")
