import math

import fasta

class BiasReport (object):
  def __init__(self, candidate_fasta, reference_fasta, donor_vcf, log, buckets):
    '''
      @candidate_fasta: results from alignment
      @reference_fasta: what was aligned against
      @donor_vcf: VCF file used to generate reads
      @log
      @buckets
      TODO only handles SNPs
    '''
    reference = fasta.Fasta( fasta.FastaReader( reference_fasta ) )
    pos = 0
    self.reference_histogram = [ 0 ] * buckets
    self.error_histogram = [ 0 ] * buckets
    self.stats = { 'unmapped': 0, 'mapped': 0 }
    while pos < candidate_fasta.length:
      reference_base = reference.base_at( pos )
      if reference_base is None:
        break # reference is done
      # check for SNP here TODO any variation
      if pos in donor_vcf.snp_map:
        counts = candidate_fasta.count( pos )
        total = sum( [ counts[base] for base in counts ] )
        if total > 0:
          self.stats['mapped'] += 1
          donor_base = donor_vcf.snp_list[donor_vcf.snp_map[pos]]['alt']
          reference_percent = 1. * counts[reference_base] / total
          error_percent = 1. * ( total - counts[reference_base] - counts[donor_base] ) / total
          self.add_to_stats( reference_percent, error_percent, buckets )
          #print "pos", pos, "counts", counts, "ref", reference_base, "vcf", donor_vcf.snp_list[donor_vcf.snp_map[pos]], "ref%", reference_percent, "err%", error_percent
        else: # no coverage
          self.stats['unmapped'] += 1 # is this bias?
      #if candidate_base != reference_base and candidate_base != 'N':
      #  target_vcf.snp( pos, reference_base, candidate_base ) # mutation
      pos += 1

  def add_to_stats( self, reference, error, buckets ):
    # find reference bucket
    bucket_size = 1. / buckets
    reference_bucket = int( math.floor( reference / bucket_size ) )
    if reference_bucket >= len(self.reference_histogram):
      reference_bucket = len(self.reference_histogram) - 1
    error_bucket = int( math.floor( error / bucket_size ) )
    if error_bucket >= len(self.error_histogram):
      error_bucket = len(self.error_histogram) - 1
    self.reference_histogram[reference_bucket] += 1
    self.error_histogram[error_bucket] += 1
