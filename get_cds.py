# -*- coding: utf-8 -*-
"""
Created on Mon Aug  6 09:57:06 2018

@author: yaoxinzhi
"""

import re
import os
import argparse

def gtf_cds(gtf, positon_file):
    print ('Start processing gtf files')
    p = re.compile('.+?"(.+?)"')
    
    result = open(positon_file, 'w')
#    result.write('#{0}\t{1}\t{2}\t{3}\t{4}\n'.format('chromosome', 'gene_id', 'began', 'end', 'strand'))
    with open(gtf) as f:
        chr_num = 0
        for line in f:
            if not str.startswith(line, '#'):
                                  l = line.split('\t')
                                  if  l[0] != chr_num:
                                      chr_num = l[0]
                                      result.write('#{0}\tdna:chromosome\n'.format(chr_num))
                                  if l[2] == 'CDS':
                                     
                                      gene_id = p.findall(l[-1].split(';')[0])[0]
                                      result.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format(l[0], gene_id, l[3], l[4], l[6]))
    result.close()
    print ('CDS location file: "{0}" Get!!'.format(positon_file))  
    
def read_fa(fa_file):
    with open(fa_file) as f:
        f.readline()
        seq = f.read().replace('\n', '')
    return seq
    
def split_fa(fa_file, prefix, out_dic):
    
    print ('Start splitting the reference genome file')
    if not os.path.exists(out_dic):
        os.mkdir(out_dic)
    with open(fa_file) as f:
        first_line = True
        for line in f:
            if str.startswith(line, '>'):
                chr_num = line.split()[0][1:]
                print ('Processing sequence of No.{0} chromosome'.format(chr_num))
#                out = out_dic + '/' + prefix + '_' + chr_num + '.fa'
                out = '{0}/{1}_{2}.fa'.format(out_dic, prefix, chr_num)
                if first_line:
                    out_file = open(out, 'w')
                    out_file.write(line)
                if not first_line:
                    out_file.close()
                    out_file = open(out, 'w')
                    out_file.write(line)
            else:
                out_file.write(line)
        out_file.close()
    print ('Reference genome file segmentation completed\n')
     
def gtf2cds(fa_prefix, cds_index_file, out_file):
    print ('Start extracting cds sequences')
    with open(out_file, 'w') as f1:
        with open(cds_index_file) as f:
            gene_dic = {}
            for line in f:
                l = line.replace('\n', '').split()
                if str.startswith(line, '#'):
                                  chr_num = l[0][1:]
                                  print ('Processing No.{0} chromosome'.format(chr_num))
                                  fa_file = '{0}_{1}.fa'.format(fa_prefix, chr_num)
                                  fa_seq = read_fa(fa_file)
                                  gene_dic = {}
                else:
                    
                    gene_id = l[1]
                    if gene_id in gene_dic.keys():
                        gene_dic[gene_id] += 1
                    else:
                        gene_dic[gene_id] = 1
                    if gene_dic[gene_id] > 1:
                        gene_id = '{0}.{1}'.format(gene_id, gene_dic[gene_id])
                    
                    began = int(l[2])
                    end = int(l[3])
                    f1.write('{0}\t{1}\t{2}\t{3}'.format(l[0], gene_id, l[2], l[3]))
                    f1.write('\n')
                    f1.write(fa_seq[began+1:end+1] + '\n')

def main():
    
    parse = argparse.ArgumentParser()
    parse.add_argument('-g', '--gtf', help='Specify the path of GTF file',
                       dest='gtf', required=True)
    parse.add_argument('-f', '--fna', help='Specify the path of Reference FA file',
                       dest='fna', required=True)
    parse.add_argument('-c', '--cds_file', help='Specify the cds index output file, default "cds_index.csv"',
                       dest='cds_index_file', default='cds_index.csv')
    parse.add_argument('-s', '--fa_split', help='Specify the fa split file output dic, default"Split"',
                       dest='fa_output', default='Split')
    parse.add_argument('-p', '--fa_prefix', help='Specify fa split file prefix, default "Homo_sapiens"',
                       dest='fa_presix', default='Homo_sapiens')
    parse.add_argument('-o', '--output', help='Specify cds extraction result output file, default "result.txt"',
                       dest='result_output', default='result.txt')
    args = parse.parse_args()
#    
#    gtf = 'D:/yaoxinzhi/CDS/data/Homo_sapiens.GRCh38.89.chr.gtf'
#    fna = 'D:/yaoxinzhi/CDS/data/Homo_sapiens.GRCh38.dna.toplevel.fa'
#    cds_indes_file = 'cds_index.csv'
#    fa_output = 'test'
#    fa_prefix = 'Homo_sapiens'
#    result_output = 'result.fa'
#    
    gtf_cds(args.gtf, args.cds_indes_file)
    
    split_fa(args.fna, args.fa_prefix, args.fa_output)
    
    fa_split_prefix = '{0}/{1}'.format(args.fa_output, args.fa_prefix)
    gtf2cds(fa_split_prefix, args.cds_indes_file, args.result_output)

    
if __name__ == '__main__':
    main()