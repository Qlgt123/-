#!/usr/bin/env python3：声明了脚本应该使用Python 3来解释执行。
from pyfaidx import Fasta#从 pyfaidx 模块中导入了 Fasta 类，用于读取FASTA格式的文件。
import sys#导入了 sys 模块，用于读取命令行参数。

if len(sys.argv) != 3:#这段代码检查了命令行参数的数量是否为3（包括脚本本身的名称）。如果不等于3，则会打印出正确的用法说明，并退出程序。这里sys.argv是一个列表，包含了命令行传递给Python脚本的所有参数。
    print("用法：python extract_HGT_sequences.py verified_HGT.id.txt query_oneline.fasta")
    sys.exit(1)
id_file = sys.argv[1]
fasta_file = sys.argv[2]
#id_file 和 fasta_file 分别从命令行参数中获取ID列表文件和基因组FASTA文件的路径。
output_file = "verified_HGT_sequences.fasta"#output_file 定义了输出文件的名称，即提取出的序列将写入这个文件中。
# 读取基因组
genome = Fasta(fasta_file)#使用 Fasta 类读取FASTA格式的基因组文件，genome 变量包含了基因组的序列信息。
# 读取ID列表
with open(id_file) as f:
    ids = [line.strip() for line in f if line.strip()]#打开ID列表文件，逐行读取其内容。line.strip() 用于去除每行首尾的空白字符，包括换行符。
#列表推导式 [line.strip() for line in f if line.strip()] 只保留了非空行的内容，存储在 ids 列表中。
# 提取序列
with open(output_file, 'w') as out:#打开输出文件，并将提取出的序列写入其中。
    for name in ids:#对 ids 列表中的每一个ID进行处理
        try:
            parts = name.split('---')#parts = name.split('---') 将ID字符串根据’—'分割成两个部分，第一个部分是contig名称，第二个部分是序列的起始和结束位置。
            contig = parts[0]# 获取contig名称
            start, end = map(int, parts[1].split('_')[:2])#start, end = map(int, parts[1].split('_')[:2]) 从ID的第二部分提取起始和结束位置，并将其转换为整数类型点击链接查看和 Kimi 的对话 https://www.kimi.com/share/19b01bcd-17d2-80b3-8000-00004dbb1082
            seq = genome[contig][start:end].seq#seq = genome[contig][start:end].seq 使用这些位置信息从基因组中提取对应的序列。
            out.write(f'>{name}\n{seq}\n')#out.write(f'>{name}\n{seq}\n') 将提取出的序列以FASTA格式写入输出文件中。点击链接查看和 Kimi 的对话 https://www.kimi.com/share/19b01bf3-2492-8488-8000-0000da2abf57
        except Exception as e:
            print(f"[跳过] {name} 出错：{e}")

print(f"[完成] 已提取 {len(ids)} 条序列到 {output_file}")
