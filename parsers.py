import re
from collections import OrderedDict

class Variant:
    def __init__(self, line, samples):
        columns = line.strip().split('\t')
        self.CHROM = columns[0]
        self.POS = int(columns[1])
        self.ID = columns[2]
        self.REF = columns[3]
        self.ALT = columns[4].split(',') if columns[4] != '.' else []
        self.QUAL = columns[5] if columns[5] != '.' else None
        self.FILTER = columns[6].split(';') if columns[6] != '.' else []
        
        self.INFO = {}
        info_str = columns[7]
        if info_str != '.':
            for ent in info_str.split(';'):
                if '=' in ent:
                    k, v = ent.split('=', 1)
                    self.INFO[k] = v
                else:
                    self.INFO[ent] = True
        
        self.FORMAT = []
        self.samples = OrderedDict()
        
        if len(columns) > 8:
            self.FORMAT = columns[8].split(':')
            for i, sample in enumerate(samples):
                sample_data = columns[9 + i] if 9 + i < len(columns) else '.'
                sample_values = sample_data.split(':')
                sample_dict = {}
                for j, key in enumerate(self.FORMAT):
                    sample_dict[key] = sample_values[j] if j < len(sample_values) else '.'
                self.samples[sample] = sample_dict
    
    def genotype(self, sample_name):
        return self.samples.get(sample_name, {})

    def genotypes(self):
        return {sample: info.get('GT', '.') for sample, info in self.samples.items()}


class VCFParser:
    def __init__(self, vcf_file):
        self.vcf_file = vcf_file
        
        self.ref = ''
        self.fformat = ''
        self.contig = OrderedDict()
        self.samples = []
        
        self.INFO = OrderedDict()
        self.FORMAT = OrderedDict()
        self.FILTER = OrderedDict()
        self.ALT = OrderedDict()
        self.other = []
        
        self._parse_header()
        
    def _parse_header(self):

        info_pattern = re.compile(
            r'##INFO=<ID=(?P<id>[^,>]+),'
            r'Number=(?P<number>[^,>]+),'
            r'Type=(?P<type>[^,>]+),'
            r'Description="(?P<desc>(?:[^"]|\\")*)"'
            r'(?:,.*)?'
            r'>'
        )
        
        format_pattern = re.compile(
            r'##FORMAT=<ID=(?P<id>[^,>]+),'
            r'Number=(?P<number>[^,>]+),'
            r'Type=(?P<type>[^,>]+),'
            r'Description="(?P<desc>(?:[^"]|\\")*)"'
            r'(?:,.*)?'
            r'>'
        )
        
        filter_pattern = re.compile(
            r'##FILTER=<ID=(?P<id>[^,>]+),'
            r'Description="(?P<desc>(?:[^"]|\\")*)"'
            r'(?:,.*)?'
            r'>'
        )
        
        alt_pattern = re.compile(
            r'##ALT=<ID=(?P<id>[^,>]+),'
            r'Description="(?P<desc>(?:[^"]|\\")*)"'
            r'(?:,.*)?'
            r'>'
        )

    
        contig_pattern = re.compile(
            r'##contig=<'
            r'[^>]*\bID=(?P<chrom>[^,>]+)'                 
            r'(?:[^>]*\blength=(?P<length>\d+))?'          
            r'(?:[^>]*\bassembly=(?P<assembly>[^,>]+))?'     
            r'[^>]*>'
        )
    
        fileformat_pattern = re.compile(r'^##fileformat=')
        reference_pattern = re.compile(r'^##reference=')
    
        info_idx = format_idx = filter_idx = alt_idx = contig_idx = 0
        
        with open(self.vcf_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                if line.startswith('##'):
                    if fileformat_pattern.match(line):
                        self.fformat = line
    
                    elif reference_pattern.match(line):
                        self.ref = line
    
                    elif info_pattern.match(line):
                        m = info_pattern.match(line)
                        self.INFO[m.group('id')] = {
                            'number': m.group('number'),
                            'type': m.group('type'),
                            'description': m.group('desc'),
                        }
                        info_idx += 1
    
                    elif format_pattern.match(line):
                        m = format_pattern.match(line)
                        self.FORMAT[m.group('id')] = {
                            'number': m.group('number'),
                            'type': m.group('type'),
                            'description': m.group('desc'),
                        }
                        format_idx += 1
    
                    elif filter_pattern.match(line):
                        m = filter_pattern.match(line)
                        self.FILTER[m.group('id')] = {
                            'description': m.group('desc'),
                        }
                        filter_idx += 1
    
                    elif alt_pattern.match(line):
                        m = alt_pattern.match(line)
                        self.ALT[m.group('id')] = {
                            'description': m.group('desc'),
                        }
                        alt_idx += 1
    
                    elif contig_pattern.match(line):
                        m = contig_pattern.match(line)
                        self.contig[m.group('chrom')] = {
                            'length': int(m.group('length')) if m.group('length') else None,
                            'assembly': m.group('assembly') if m.group('assembly') else None,
                        }
                        contig_idx += 1
    
                    else:
                        self.other.append(line)
    
                elif line.startswith('#CHROM'):
                    headers = line.split('\t')
                    self.samples = headers[9:]
                    break  #метаинформация закончилась, дошли до заголовка, дальше тело, надо переходить к чтению вариантов
    
    def parse_variants(self):
        with open(self.vcf_file, 'r') as f:
            in_body = False
            for line in f:
                if line.startswith('#CHROM'):
                    in_body = True
                    continue
                if not in_body:
                    continue
                line = line.strip()
                if not line:
                    continue
                yield Variant(line, self.samples)
    
    def get_variants(self, chrom, start, length):
        end = start + length
        for variant in self.parse_variants():
            if variant.CHROM == chrom and start <= variant.POS < end:
                yield variant


class FASTAParser:
    def __init__(self):
        self.sequences = {}

    def parse(self, ref_genome_path):
        current_chr = None
        buffer = []
        with open(ref_genome_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if current_chr is not None:
                        self.sequences[current_chr] = ''.join(buffer).upper()
                    current_chr = line[1:].split()[0]
                    buffer = []
                else:
                    if current_chr is not None:
                        buffer.append(line)
            if current_chr is not None:
                self.sequences[current_chr] = ''.join(buffer).upper()

    def get_sequence(self, chrom, start, end):
        """Возвращает последовательность из региона указанной длины"""
        if chrom not in self.sequences:
            return None
        seq = self.sequences[chrom]
        if not (1 <= start <= end <= len(seq)):
            raise ValueError(f"Проверьте координаты: start={start}, end={end}, длина={len(seq)}")
        return seq[start-1:end]

    def chr_names(self):
        return list(self.sequences.keys())
