#!/usr/bin/env python3
import argparse
import sys
import os

from parsers import VCFParser, FASTAParser
from samplers import RandomSampler, ConsGenerator

def main():
    parser = argparse.ArgumentParser(
        description='Генератор последовательностей из VCF и FASTA файла'
    )
    
    parser.add_argument('--vcf', required=True, help='Путь к VCF файлу')
    parser.add_argument('--fasta', required=True, help='Путь к референсному FASTA файлу')
    parser.add_argument('-O', required=True, help='Выходной файл')
    parser.add_argument('-L', required=True, type=int, metavar='INT', help='Длина последовательности')
    parser.add_argument('-c', required=True, type=int, metavar='INT', help='Количество последовательностей')
    parser.add_argument('--minAF', required=True, type=float, metavar='FLOAT', help='Минимальная частота по AF')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.vcf):
        print(f"VCF не найден в {args.vcf}", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(args.fasta):
        print(f"FASTA не найден в {args.fasta}", file=sys.stderr)
        sys.exit(1)
    
    try:
        vcf_parser = VCFParser(args.vcf)
        fasta_parser = FASTAParser()
        fasta_parser.parse(args.fasta)
        
        chr_name = fasta_parser.chr_names()[0]
        
        random_sampler = RandomSampler(fasta_parser)
        cons_generator = ConsGenerator(chr_name=chr_name, vcf_parser=vcf_parser,
                                       fasta_parser=fasta_parser, random_sampler=random_sampler,
                                       num_seqs=args.count, seq_length=args.length,min_af=args.minAF)
        
        sample_variants_dict, sequences = cons_generator.sample_variants()
        
        results = []
        for sample, variants_list in sample_variants_dict.items():
            for start_pos, _ in sorted(sequences.items()):
                consensus_seq = cons_generator.generate_consensus(
                    sample, chr_name, start_pos, args.length, variants_list
                )
                header = f">{sample}_consensus_{start_pos}_{start_pos + args.length - 1}"
                results.append(header)
                results.append(consensus_seq)
        
        with open(args.output, 'w') as f:
            for line in results:
                f.write(line + '\n')
        
    except Exception as e:
        print(f"{e}", file=sys.stderr)
        sys.exit(1)
    
    sys.exit(0)

if __name__ == '__main__':
    main()
