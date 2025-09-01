import random
import re

class RandomSampler:
    def __init__(self, fasta_parser, random_seed=None):
        self.fasta_parser = fasta_parser
        self.random_seed = random_seed
        
    def generate_seqs(self, chrom, num_sequences, sequence_length):
        """
        Генерирует случайные последовательности определенной длины.
        Возвращает словарь {start: sequence}
        """
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
    def _is_valid_alt(alt_list):
        pattern = re.compile(r'^[ATCG]+$')
        return all(pattern.match(alt) for alt in alt_list)

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
        recording = False
        for variant in self.vcf_parser.parse_variants():
            if variant.CHROM != self.chr_name:
                continue
            if variant.POS >= min_start:
                recording = True
            if not recording:
                continue
            if variant.POS > max_start + self.seq_length:
                break
            if not self._is_valid_alt(variant.ALT):
                continue
            for sample in self.vcf_parser.samples:
                gt = variant.genotype(sample).get('GT', './.')
                if gt not in ['0/0', '0|0', './.']:
                    samples_dict[sample].append(variant)
        return samples_dict, seqs

    def generate_consensus(self, sample, chrom, start, length, sample_variants):
        ref_seq = self.fasta_parser.get_sequence(chrom, start, start + length - 1)
        variants_in_region = [v for v in sample_variants if start <= v.POS < start + length]
        consensus_seq = self.apply_variants(ref_seq, variants_in_region, sample)
        return consensus_seq

    def parse_genotype(self, gt_string):
        if gt_string in ['./.', '.|.', '.']:
            return None, None

        if '|' in gt_string:
            alleles = gt_string.split('|')
            phased = True
        else:
            alleles = gt_string.split('/')
            phased = False

        return [int(a) if a.isdigit() else None for a in alleles], phased

    def apply_variants(self, ref_seq, variants, sample):
        sequence = list(ref_seq)
        cumulative_offset = 0

        for variant in sorted(variants, key=lambda v: v.POS):
            adjusted_pos = variant.POS - 1 + cumulative_offset
            if adjusted_pos < 0 or adjusted_pos >= len(sequence):
                continue

            gt_str = variant.genotype(sample).get('GT', './.')
            alleles, _ = self.parse_genotype(gt_str)
            if not alleles or all(a == 0 or a is None for a in alleles):
                continue

            alt_idx = None
            alt_allele = None

            af_list = [
                float(x) for x in variant.INFO.get('AF', '').split(',') if x
            ] if 'AF' in variant.INFO else None

            nonref_alleles = [a for a in alleles if a and a > 0]
            if not nonref_alleles:
                continue

            if len(set(nonref_alleles)) == 1:
                candidate_idx = nonref_alleles[0] - 1
                if af_list and candidate_idx < len(af_list):
                    if af_list[candidate_idx] >= self.min_af:
                        alt_idx = candidate_idx
                else:
                    alt_idx = candidate_idx
            else:
                if af_list:
                    valid_alleles = [a for a in nonref_alleles if 0 < (a - 1) < len(af_list) and af_list[a-1] >= self.min_af]
                    if valid_alleles:
                        best = max(valid_alleles, key=lambda idx: af_list[idx - 1])
                        alt_idx = best - 1
                else:
                    alt_idx = nonref_alleles[0] - 1

            if alt_idx is not None and 0 <= alt_idx < len(variant.ALT):
                alt_allele = variant.ALT[alt_idx]
                ref_len = len(variant.REF)
                alt_len = len(alt_allele)
                sequence[adjusted_pos:adjusted_pos + ref_len] = list(alt_allele)
                cumulative_offset += alt_len - ref_len

        return ''.join(sequence)
