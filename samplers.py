import random
import re

class RandomSampler:
    def __init__(self, fasta_parser, random_seed=None):
        self.fasta_parser = fasta_parser
        self.random_seed = random_seed

    def generate_seqs(self, chrom, num_sequences, sequence_length):
        if chrom not in self.fasta_parser.chr_names():
            raise ValueError(f"{chrom} not in fasta_parser.chr_names()")

        chr_length = len(self.fasta_parser.sequences[chrom])
        if sequence_length > chr_length:
            raise ValueError(
                f"Specified length {sequence_length} "
                f"can't be put into {chrom} with length {chr_length}"
            )

        max_right = chr_length - sequence_length + 1
        if num_sequences > max_right:
            raise ValueError(
                "More sequences than allowed to avoid non unique start positions"
            )

        if self.random_seed is not None:
            random.seed(self.random_seed)

        starts = random.sample(range(1, max_right + 1), k=num_sequences)
        sequences = {
            s: self.fasta_parser.get_sequence(chrom, s, s + sequence_length - 1)
            for s in starts
        }
        return sequences


class ConsGenerator:
    def __init__(self, chr_name, vcf_parser, fasta_parser, random_sampler, num_seqs, seq_length, min_af=0.0):
        self.chr_name = chr_name
        self.vcf_parser = vcf_parser
        self.fasta_parser = fasta_parser
        self.random_sampler = random_sampler
        self.num_seqs = num_seqs
        self.seq_length = seq_length
        self.min_af = min_af

    @staticmethod
    def is_valid_alt(alt_list):
        pattern = re.compile(r'^[ATCG]+$')
        return all(pattern.match(alt) for alt in alt_list)

    def parse_genotype(self, genotype):
        if genotype in ['./.', '.|.', '.']:
            return None, None
        
        if '|' in genotype:
            alleles = genotype.split('|')
            phased = True
        else:
            alleles = genotype.split('/')
            phased = False
        
        try:
            allele_indices = [int(a) for a in alleles if a.isdigit()]
        except ValueError:
            return None, None
        
        return allele_indices, phased

    def sample_variants(self):
        samples_dict = {sample: [] for sample in self.vcf_parser.samples}
        seqs = self.random_sampler.generate_seqs(
            chrom=self.chr_name,
            num_sequences=self.num_seqs,
            sequence_length=self.seq_length
        )

        starts = sorted(seqs.keys())
        min_start = starts[0]
        max_start = starts[-1]

        for variant in self.vcf_parser.parse_variants():
            if variant.CHROM != self.chr_name:
                continue
            
            if not (min_start <= variant.POS <= max_start + self.seq_length):
                continue

            if not self.is_valid_alt(variant.ALT):
                continue

            if not self.passes_af_filter(variant):
                continue

            for sample in self.vcf_parser.samples:
                gt = variant.genotype(sample).get('GT', './.')
                alleles, _ = self.parse_genotype(gt)
                
                if alleles and any(a > 0 for a in alleles):
                    samples_dict[sample].append(variant)

        return samples_dict, seqs

    def passes_af_filter(self, variant):
        if 'AF' not in variant.INFO:
            return True
        
        af_values = variant.INFO['AF'].split(',')
        try:
            return any(float(af) >= self.min_af for af in af_values if af)
        except ValueError:
            return True

    def generate_consensus(self, sample, chrom, start, length, sample_variants):
        ref_seq = self.fasta_parser.get_sequence(chrom, start, start + length - 1)
        
        variants_in_region = [
            v for v in sample_variants 
            if start <= v.POS < start + length
        ]
        
        if not variants_in_region:
            return ref_seq
        
        return self.apply_variants(ref_seq, variants_in_region, sample, start)

    def choose_best_allele(self, variant, genotype):
        alleles, phased = self.parse_genotype(genotype)
        if not alleles:
            return None
        
        alt_alleles = [a - 1 for a in alleles if a > 0]
        if not alt_alleles:
            return None
        
        valid_alleles = []
        af_values = []
        
        if 'AF' in variant.INFO:
            try:
                af_values = [float(x) for x in variant.INFO['AF'].split(',') if x]
            except ValueError:
                af_values = []
        
        for allele_idx in alt_alleles:
            if allele_idx < len(variant.ALT):
                if af_values and allele_idx < len(af_values):
                    if af_values[allele_idx] >= self.min_af:
                        valid_alleles.append((allele_idx, af_values[allele_idx]))
                else:
                    valid_alleles.append((allele_idx, 1.0))
        
        if not valid_alleles:
            return None
        
        if len(set(alleles)) == 1:
            return valid_alleles[0][0]
        
        best_allele = max(valid_alleles, key=lambda x: x[1])
        return best_allele[0]

    def apply_variants(self, ref_seq, variants, sample, region_start):
        sequence = list(ref_seq)
        
        sorted_variants = sorted(variants, key=lambda v: v.POS, reverse=True)
        
        for variant in sorted_variants:
            rel_pos = variant.POS - region_start
            
            if rel_pos < 0 or rel_pos >= len(sequence):
                continue
            
            gt = variant.genotype(sample).get('GT', './.')
            chosen_allele_idx = self.choose_best_allele(variant, gt)
            
            if chosen_allele_idx is not None:
                ref_len = len(variant.REF)
                alt_seq = variant.ALT[chosen_allele_idx]
                
                sequence[rel_pos:rel_pos + ref_len] = list(alt_seq)
        
        return ''.join(sequence)