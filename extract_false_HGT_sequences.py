#!/usr/bin/env python3这行代码指定了使用 Python 3 解释器来运行脚本。env 命令会根据系统环境找到合适的解释器。
"""
extract_false_HGT_sequences.py
仅依赖 pyfaidx，用于从核基因组中提取
被误判为“水平基因转移”的染色体片段（False HGT）
用法：
python extract_false_HGT_sequences.py verified_false_HGT.id.txt query_oneline.fasta
"""
import sys
from pyfaidx import Fasta#sys 库用于与 Python 解释器进行交互，特别是处理命令行参数。pyfaidx.Fasta 是一个用于读取 FASTA 格式文件的类，提供了方便的接口来访问序列数据。
if len(sys.argv) != 3:
    print("用法：python extract_false_HGT_sequences.py verified_false_HGT.id.txt query_oneline.fasta")
    sys.exit(1)
id_file   = sys.argv[1]          # verified_false_HGT.id.txt
fasta_file= sys.argv[2]          # query_oneline.fasta（核基因组）
out_file  = "verified_false_HGT_sequences.fasta"
genome = Fasta(fasta_file)#这行代码使用 pyfaidx.Fasta 类创建了一个 Fasta 对象，该对象表示核基因组序列文件。通过这个对象，可以方便地访问和操作文件中的序列数据。
with open(id_file) as f:#使用 with open 语句打开 id_file 文件。
    ids = [line.strip() for line in f if line.strip()]
    #通过列表推导式读取文件中的每一行，并使用 strip() 方法去除每行前后的空白字符（包括换行符）。只保留非空行，存储在 ids 列表中。
with open(out_file, 'w') as out:#使用 with open 语句打开 out_file 文件，准备写入数据。
    for name in ids:
        try:
            parts   = name.split('---')#使用 split('---') 方法将 ID 分割为两部分：contig 和后续的部分。
            contig  = parts[0]
            start, end = map(int, parts[1].split('_')[:2])#使用 split('_')[:2] 方法将后续部分分割为起始位置和结束位置。使用 map(int, ...) 方法将起始位置和结束位置转换为整数
            seq     = genome[contig][start:end].seq#通过 genome[contig][start:end].seq 提取指定染色体片段的序列
            out.write(f'>{name}\n{seq}\n')#使用 out.write 方法将序列的 FASTA 格式写入输出文件
        except Exception as e:
            print(f"[跳过] {name} 出错：{e}")
print(f"[完成] 已提取 {len(ids)} 条误判HGT序列到 {out_file}")
#该脚本的主要功能是从一个核基因组 FASTA 文件（query_oneline.fasta）中提取被误判为水平基因转移（HGT）的染色体片段。需要提供一个包含这些片段 ID 的文件（verified_false_HGT.id.txt），脚本会根据这些 ID 提取对应的序列，并将提取的序列写入一个新的 FASTA 文件（verified_false_HGT_sequences.fasta）。每个序列 ID 应该格式为 contig---start_end，其中 contig 是染色体的名称，start 和 end 是需要提取的片段的起始和结束位置。
