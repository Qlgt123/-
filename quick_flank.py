# 快速提取HGT上下游10kb
# 输入：merged.bed 和 query_oneline.fasta
# 输出：HGT_10kb_flank.fasta + HGT_10kb_summary.tsv

import re
from Bio import SeqIO

genome = SeqIO.to_dict(SeqIO.parse("query_oneline.fasta", "fasta"))
flank = 10000

with open("merged.bed") as bed, \
     open("HGT_10kb_flank.fasta", "w") as fasta_out, \
     open("HGT_10kb_summary.tsv", "w") as tsv_out:

    tsv_out.write("contig\tHGT_start\tHGT_end\tupstream_start\tupstream_end\tdownstream_start\tdownstream_end\tup_GC\tdown_GC\n")

    for line in bed:
        c, s, e = line.strip().split("\t")[:3]
        s, e = int(s), int(e)
        up_start = max(0, s - flank)
        down_end = e + flank

        up_seq   = str(genome[c].seq[up_start:s]).upper()
        down_seq = str(genome[c].seq[e:down_end]).upper()

        def gc(seq): return round((seq.count("G") + seq.count("C")) / len(seq) * 100, 2) if seq else 0

        fasta_out.write(f">{c}:{s}-{e}_upstream\n{up_seq}\n")
        fasta_out.write(f">{c}:{s}-{e}_downstream\n{down_seq}\n")

        tsv_out.write("\t".join(map(str, [c, s, e, up_start, s, e, down_end, gc(up_seq), gc(down_seq)])) + "\n")

print("✅ 完成：HGT_10kb_flank.fasta 和 HGT_10kb_summary.tsv 已生成")

