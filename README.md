1、extract_HGT_sequences.py提取正确的水平基因转移
2、extract_false_HGT_sequences.py提取错误的水平基因转移
3、quick.flank.py水平基因转移的位置（基因组哪些位置）的上游与下游（构成水平基因转移的上下游序列的特点），是不是有些位置容易发生水平基因转移，前面是哪些基因发生了水平基因转移，这里是哪些位置发生了水平基因转移，提取出序列
4、quick_flank_false_HGT与3类似，但是是错误的水平基因转移位置
5、3、4的代码要改位置：（存放数据的位置）
GENOME_FA   = "/home/huajin/newdisk_1/sxk/songwei/organelle_genome_contamination_verification/shuidao/10000/fchlo/query_oneline.fasta"
CANDIDATE_FA= "/home/huajin/newdisk_1/sxk/songwei/organelle_genome_contamination_verification/shuidao/10000/fchlo/extracted_HTG_updown_region.fasta"
