import os
import sys
from numpy import *
import math
import numpy as np
import multiprocessing
from numpy import *
import math
import time
import datetime
from datetime import datetime, date
from numpy import *
import math
from PyPDF2 import PdfFileMerger
import numpy as np


# (base) root@dell-server:/home/newdisk_dell_6/wilson/organelle_genome_contamination_verification/human# python minimap2_result_visulation-main-verify-with-TGS.py   human_mito.fasta  hg38-rm-mt.fa   human_mito_identify_result_3  1000  1000 SRR29483267_subreads.fasta


#(base) root@dell-server:/home/newdisk_dell_6/wilson/organelle_genome_contamination_verification/yumi# docker run -v /var/run/docker.sock:/var/run/docker.sock -v `pwd`:/data -w /data   songweidocker/chlomito:v1  python minimap2_result_visulation-main-verify-with-TGS-noplot.py  yumi_chlo.fasta   yumi_genome.fasta   yumi_chlo_identify_result_2  100  1000   SRR15447414_subreads.fasta

database = sys.argv[1]
query = sys.argv[2]
output = sys.argv[3]

HGT_length_cutoff = sys.argv[4]
HGT_up_down_length = sys.argv[5]

TGS_file = sys.argv[6]

HGT_distence = sys.argv[7]




os.system('rm -rf %s && mkdir %s'%(output,output))
os.chdir('%s'%(output))

os.system('cp ../%s ../%s ./'%(database,query))
os.system('cut -d " " -f 1 %s > database_oneline.fasta'%(database))
os.system('cut -d " " -f 1 %s > query_oneline.fasta'%(query))

database = "database_oneline.fasta"
query = "query_oneline.fasta"


os.system('grep  ">" %s | cut -d ">" -f 2 > contig.id '%(query))
app_tmp_0 = open("./contig.id","r").read().strip().split("\n")



os.system('samtools faidx %s'%(database))
contig_length_dic = {}
for eles in open("database_oneline.fasta.fai").read().strip().split("\n") :
    ele = eles.strip().split("\t")
    contig_length_dic[ ele[0].strip() ] = ele[1].strip()  
print(contig_length_dic)


os.system(' minimap2  %s  %s  -t 60 > mito-genome.blast.result.paf  '%(database,query))
os.system(' minimap2  %s  %s  -t 60 | awk  \'{if ($10/$11 > 0.5) print }\' - |  awk \'{print $1"\t"$2"\t"$6"\t"$7"\t""100""\t"$10"\t"$3"\t"$4"\t"$8"\t"$9}\' -  > mito-genome.blast.result.txt  '%(database,query))

filein = "mito-genome.blast.result.txt" # blastn result

def each_contig_in_database(target_contig) :   
    os.system('rm -rf %s && mkdir %s'%(target_contig,target_contig))
    os.chdir('%s'%(target_contig)) 
    os.system('cp ../%s ./'%(filein) )
    filein2 = filein + "2"
    os.system(' sort -k 1  -n -k 7 -n -k 8 %s  |awk \'{if ($6 > 10 && $3 == "%s"  )  print }\'  > %s  '%(filein,target_contig, filein2))
    os.system(' awk  -F "\t"  \'{print $0 > $1".txt"}\'  %s   '%(filein2))

    list0 = os.popen('ls ').read().strip().split("\n")
    list0.remove(filein)
    list0.remove(filein2)
    print(list0)
    
    def mapped_length_calculate(tmpfile) :
        eles = open(tmpfile).read().strip().split("\n")
        bb = []
        i = 0
        for ele in eles :
            start = ele.strip().split("\t")[8]
            end = ele.strip().split("\t")[9]
            aa = list( np.arange(int(start),int(end) ))
            bb.extend(aa)
            i = i + 1
        cc_tmp = list(set(bb))
        cc = sorted(cc_tmp) 
        for ele2 in cc :
            with open(tmpfile + "-tmp", "a+") as f :
                f.write(  tmpfile.strip().split(".txt")[0].strip() + "\t" + str(ele2) +  "\n" )
                
    for ele3 in list0 :
        mapped_length_calculate(ele3)
    os.system('cat *-tmp  | sed  \'1i seq\tposition \' > final_result-plot.txt')
    os.system('cp ../../minimap2_result_visulation.R  ./')
    input_1 = "final_result-plot.txt"
    contig_name = target_contig
    contig_length_input = contig_length_dic[target_contig]
    output_1 = target_contig + "-1.pdf"
    output_2 = target_contig + "-2.pdf"
    #os.system(' docker run  -v /var/run/docker.sock:/var/run/docker.sock -v `pwd`:/data -w /data  minimap2_result_visualization_docker_image  Rscript minimap2_result_visulation.R  %s %s %s %s %s '%( input_1, contig_name, contig_length_input, output_1, output_2  ))

    os.chdir('../')

os.system('grep  ">" %s | cut -d ">" -f 2 > database_contig.id '%(database))
database_tmp_0 = open("./database_contig.id","r").read().strip().split("\n")

for eles in database_tmp_0 :
    each_contig_in_database(eles) 





#########################################################
os.system('rm -rf TGS_verify_result && mkdir TGS_verify_result ')
os.chdir('TGS_verify_result')
os.system('cp ../mito-genome.blast.result.txt  ./')






os.system('cut -f 1,7,8 mito-genome.blast.result.txt |sort | uniq | awk \'{if (($3-$2)> %s) print  } \' |  sort -k 1,1 -k 2,2n > tmp_sorted.bed '%(HGT_length_cutoff))
os.system('bedtools merge -i tmp_sorted.bed -d %s -c 2,3 -o collapse > merged.bed '%(HGT_distence))
bed_output_1 =  HGT_length_cutoff + "_bp_HGT_"  +  HGT_up_down_length + "_bp_updown_region.bed"
os.system('awk \'{print $1"\t"$2-%s"\t"$3+%s}\'  merged.bed  >  %s  '%(HGT_up_down_length,HGT_up_down_length,bed_output_1))



genome_dic = {}
sequences = open("../query_oneline.fasta").read().strip().split(">")[1:]
for eles in sequences :
    ele = eles.strip().split("\n")
    genome_dic[ele[0]] = ele[1] 


for eles in open(bed_output_1).read().strip().split("\n"):
    ele = eles.strip().split("\t") 
    seq_len = int(ele[2]) - int(ele[1])
    seq_id = ele[0] + "---" + ele[1] + "_"  + ele[2] +  "---" +  "1" + "_" + HGT_up_down_length + "_" + str(int(HGT_up_down_length)*2) + "_" +str(int(seq_len) - int(HGT_up_down_length)*2 ) + "_" + str( int(seq_len) - int(HGT_up_down_length) ) + "_" + str(seq_len)
    sequence = genome_dic[ele[0]][int(ele[1]):int(ele[2])]
    with open("extracted_HTG_updown_region.fasta","a+") as f :
        f.write(">" +  seq_id + "\n" + sequence + "\n" ) 


os.system('ln -s ../../%s ./'%(TGS_file))



os.system(' minimap2  %s  %s  -t 60 > HGT-TGS.blast.result.paf  '%("extracted_HTG_updown_region.fasta",TGS_file))



#--------------------------------------------------------------------


os.system('awk  \'{if ($10/$11 > 0.5) print }\' HGT-TGS.blast.result.paf  > filtered_HGT-TGS.blast.result.paf  ')

os.system('rm -rf True_HGT_identify_result.tx  False_HGT_identify_result.txt  True_HGT_identify_proof.txt')
tmp = []
for eles in open("filtered_HGT-TGS.blast.result.paf").read().strip().split("\n") :
    ele_1 = eles.strip().split("\t")
    chrom_fullname = ele_1[5].strip()
    HGT_position = chrom_fullname.strip().split("---")[2].strip().split("_")
    forward_region_left = int(HGT_position[0])
    forward_region_middle = int(HGT_position[1])
    forward_region_right = int(HGT_position[2])    
    backward_region_left = int(HGT_position[3])
    backward_region_middle = int(HGT_position[4])
    backward_region_right = int(HGT_position[5])    
    if int(ele_1[7]) < int(ele_1[8]) :
        mapped_position_1 =  int(ele_1[7])
        mapped_position_2 =  int(ele_1[8])
    else :
        mapped_position_1 =  int(ele_1[8])
        mapped_position_2 =  int(ele_1[7])
    if (  mapped_position_1 < forward_region_middle - int(HGT_up_down_length) / 2  and  mapped_position_2 > forward_region_middle + int(HGT_up_down_length) / 2  )  or  ( mapped_position_1 < backward_region_middle  - int(HGT_up_down_length) / 2   and mapped_position_2 >  backward_region_middle + int(HGT_up_down_length) / 2  ) :
        with open("True_HGT_identify_result.txt","a+") as f :
            f.write(eles.strip() + "\n")
        if chrom_fullname not in tmp :
            with open("True_HGT_identify_proof.txt","a+") as f :
                f.write(eles.strip() + "\n")
            tmp.append(chrom_fullname)
    else :
        with open("False_HGT_identify_result.txt","a+") as f :
            f.write(eles.strip() + "\n")



os.system('grep ">" extracted_HTG_updown_region.fasta | cut -d ">" -f 2 > all_candidate_HGT.id.txt')
os.system('cut -f 6  True_HGT_identify_result.txt  | sort | uniq > verified_HGT.id.txt ')
os.system('grep -v -f verified_HGT.id.txt all_candidate_HGT.id.txt  > verified_false_HGT.id.txt ')














def remove_same_position_sequcnces_in_paf(paf_file, output_paf, same_sequence_keep_num) :
    os.system('sort   -k 8,8n -k 9,9n  %s >  1.paf  '%(paf_file))
    aa_1 = []
    i = 0
    for eles in open("1.paf").read().strip().split("\n") :
        ele = eles.strip().split("\t")
        # 检查是否至少有9列
        if len(ele) < 9:
            print(f"Skipping line with insufficient columns: {eles}")
            continue
        position = ele[7] + "_" + ele[8]
        if position not in aa_1   :
            with open(output_paf,"a+") as f :
                f.write(eles.strip() + "\n")
            i = i + 1
        elif  position  in aa_1  and i < int(same_sequence_keep_num) :
            i = i + 1
        else :
            i = 0
    os.system('rm -rf 1.paf')






os.system('rm -rf verified_HGT_TGS_pic &&  mkdir  verified_HGT_TGS_pic ')
os.chdir('verified_HGT_TGS_pic')
os.system('cp ../../../pafCoordsDotPlotly_modified_for_HGT_detected_in_Chrom.R ./')
os.system('cp ../verified_HGT.id.txt  ./')
for eles in open('verified_HGT.id.txt').read().strip().split("\n") :
    os.system('grep %s ../True_HGT_identify_result.txt  > tmp_0.paf '%(eles))
    remove_same_position_sequcnces_in_paf("tmp_0.paf","tmp.paf","3")
    output = eles + "_tmp_plot"
    mapped_start_position = eles.strip().split("---")[2].strip().split("_")[1]
    mapped_end_position = eles.strip().split("---")[2].strip().split("_")[4]
    os.system('Rscript ./pafCoordsDotPlotly_modified_for_HGT_detected_in_Chrom.R  -i   tmp.paf  -o  %s  -s -t -m 50 -q 50  -a %s -b %s  '%(output,mapped_start_position,mapped_end_position))
    os.system('rm -rf tmp.paf ')
files = os.listdir('./')
merger = PdfFileMerger()
for pdf in files:
    if pdf[-4:] == ".pdf":
        pdf_file="./"+pdf
        merger.append(open(pdf_file, 'rb'))
with open('verified_HGT_mapped_with_HiFi.pdf', 'wb') as fout:
    merger.write(fout)


os.system('rm -rf *_tmp_plot.pdf')
os.chdir('../')

os.system('rm -rf verified_false_HGT_TGS_pic_2 &&  mkdir  verified_false_HGT_TGS_pic_2')
os.chdir('verified_false_HGT_TGS_pic_2')
os.system('cp ../../../pafCoordsDotPlotly_modified_for_HGT_detected_in_Chrom.R ./')
os.system('cp ../verified_false_HGT.id.txt  ./')

for eles in open('verified_false_HGT.id.txt').read().strip().split("\n"):
    # 检查格式是否正确
    parts = eles.strip().split("---")
    if len(parts) < 3:
        print(f"Skipping line with insufficient parts: {eles}")
        continue

    sub_parts = parts[2].strip().split("_")
    if len(sub_parts) < 5:
        print(f"Skipping line with insufficient sub-parts: {eles}")
        continue

    mapped_start_position = sub_parts[1]
    mapped_end_position = sub_parts[4]

    os.system('grep %s ../False_HGT_identify_result.txt > tmp_0.paf ' % (eles))
    remove_same_position_sequcnces_in_paf("tmp_0.paf", "tmp.paf", "1")
    output = eles + "_tmp_plot"

    os.system(
        'Rscript ./pafCoordsDotPlotly_modified_for_HGT_detected_in_Chrom.R  -i   tmp.paf  -o  %s  -s -t -m 50 -q 50 -a %s -b %s  ' % (
        output, mapped_start_position, mapped_end_position))
    os.system('rm -rf tmp.paf')

files = os.listdir('./')
merger = PdfFileMerger()
for pdf in files:
    if pdf[-4:] == ".pdf":
        pdf_file = "./" + pdf
        merger.append(open(pdf_file, 'rb'))
with open('verified_false_HGT_mapped_with_HiFi.pdf', 'wb') as fout:
    merger.write(fout)

os.system('rm -rf *_tmp_plot.pdf')

"""
########################################################################################
#  计算每个contig的长度以及测序深度，测序深度是对这条contig的平均测序深度进行了log10处理,contigs_coverage.txt是最终文件
########################################################################################
def length_coverage_calculate(tmpfile) :
    tmpfile2 = tmpfile + "-tmp"
    os.system('  cut -d "%s" -f 3 %s >  %s  ' %("\t", tmpfile, tmpfile2 ))
    awk_command =  "'{sum+=$1}END{print sum/NR}'"
    aa = os.popen('awk %s  %s '%(awk_command,tmpfile2)).read().strip()
    bb = mean(float(aa))
    cc = os.popen('wc -l %s  '%(tmpfile2)).read().strip().split(" ")[0]
    length = str(int(cc)/1000) + "kb"
    if bb < 1 :
        genom_mean_cov = str(math.log(float(1),10))
    else :
        genom_mean_cov = str(math.log(float(bb),10))
    output = tmpfile2 + "_contigs_coverage.txt"
    with open(output,"a+")as f:
        f.write( tmpfile.strip() + "\t" + length + "\t" + str(bb)  + "\t" + str(genom_mean_cov).strip() + "\n")

def contig_depth_calculate() :
    #os.system('ln -s ../../%s  ./ &&  ln -s ../../%s  ./ '%(ngs1,ngs2))
    os.system('minimap2 -t 60 -ax  map-pb  %s    %s  | samtools view  -@ 60  -bS -F 256 - | samtools sort -@ 60  -o  sorted.aligned.bam  &&  bedtools genomecov -d  -ibam sorted.aligned.bam  > TGS_contigs.cov '%( "extracted_HTG_updown_region.fasta", TGS_file ))
    os.system('rm -rf tmp1 && mkdir tmp1 && cd tmp1 && ln -s ../TGS_contigs.cov  ./ ')
    os.chdir('tmp1')
    os.system(' awk  -F "\t"  \'{print $0 > $1".txt"}\'  %s   '%("TGS_contigs.cov"))
    list0 = os.popen('ls ').read().strip().split("\n")
    list0.remove("TGS_contigs.cov")
    ######################################################################
    pool = multiprocessing.Pool(100)
    for files in list0 :
        pool.apply_async(length_coverage_calculate,(files,))
    pool.close()
    pool.join()
    os.system('cat *_contigs_coverage.txt  >  contigs_coverage.txt  && cp contigs_coverage.txt  ../ ')
    os.chdir('../')
    

contig_depth_calculate()

os.system('echo $PWD')




#os.chdir('yumi_chlo_identify_result_2/TGS_verify_result/')

os.system('rm -rf True_HGT_pic  && mkdir True_HGT_pic ')

os.chdir('True_HGT_pic')

os.system('rm -rf *.pdf')

os.system('cp ../TGS_contigs.cov  ../verified_HGT.id.txt   ./')

# /opt/organelle_db

#os.system('cp  /opt/organelle_db/bamcov_plot_2_reference_genome.py  /opt/organelle_db/coverage_plot_2_cutoff.R ./')
os.system('cp  ../../../organelle_in_chrom_detection_plot.py ../../../organelle_in_chrom_detection_plot.R  ./')



os.system('grep -f verified_HGT.id.txt  TGS_contigs.cov  > genome.cov && python organelle_in_chrom_detection_plot.py ')




files = os.listdir('./')
merger = PdfFileMerger()
for pdf in files:
    if pdf[-4:] == ".pdf":
        pdf_file="./"+pdf
        merger.append(open(pdf_file, 'rb'))
with open('verified_HGT.pdf', 'wb') as fout:
    merger.write(fout)


os.system('rm -rf *_plot_result.pdf')
#if float(os.path.getsize("./mito.identified.plot.txt"))  > 0 :
#    os.system('mv  mito_identified_result-2.pdf  ../')
#else :
#    os.system('rm -rf mito_identified_result-2.pdf')
#        os._exit(0)

"""

