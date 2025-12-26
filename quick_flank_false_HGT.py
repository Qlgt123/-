from Bio import SeqIO

flank = 10_000

GENOME_FA   = "/home/huajin/newdisk_1/sxk/songwei/organelle_genome_contamination_verification/shuidao/10000/fchlo/query_oneline.fasta"
CANDIDATE_FA= "/home/huajin/newdisk_1/sxk/songwei/organelle_genome_contamination_verification/shuidao/10000/fchlo/extracted_HTG_updown_region.fasta"

genome = SeqIO.to_dict(SeqIO.parse(GENOME_FA, "fasta"))

false_ids = {l.strip() for l in open("verified_false_HGT.id.txt") if l.strip()}

with open("false_HGT_10kb_flank.fasta", "w") as ff, open("false_HGT_10kb_summary.tsv", "w") as tf:
    tf.write("contig\tHGT_start\tHGT_end\tupstream_start\tupstream_end\tdownstream_start\tdownstream_end\tup_GC\tdown_GC\n")
    for rec in SeqIO.parse(CANDIDATE_FA, "fasta"):
        if rec.id not in false_ids:
            continue
        c, s, e = rec.id.split("---")[0], *map(int, rec.id.split("---")[1].split("_"))
        up_start, down_end = max(0, s - flank), e + flank
        up_seq = str(genome[c].seq[up_start:s]).upper()
        down_seq = str(genome[c].seq[e:down_end]).upper()
        gc = lambda seq: round((seq.count("G") + seq.count("C")) / len(seq) * 100, 2) if seq else 0
        ff.write(f">{rec.id}_upstream\n{up_seq}\n>{rec.id}_downstream\n{down_seq}\n")
        tf.write("\t".join(map(str, [c, s, e, up_start, s, e, down_end, gc(up_seq), gc(down_seq)])) + "\n")

print("✅ false_HGT_10kb_flank.fasta & false_HGT_10kb_summary.tsv 已生成")
