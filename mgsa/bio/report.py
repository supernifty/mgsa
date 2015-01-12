import math

import fasta

class BiasReport (object):
  def __init__(self, candidate_fasta, reference_fasta, donor_vcf, log, buckets_count):
    '''
      @candidate_fasta: results from alignment
      @reference_fasta: what was aligned against
      @donor_vcf: VCF file used to generate reads
      @log
      @buckets_count
      TODO only handles SNPs
    '''
    reference = fasta.Fasta( fasta.FastaReader( reference_fasta ) )
    pos = 0
    # buckets for levels of bias in variations
    self.reference_histogram = [ 0 ] * (buckets_count+1)
    self.error_histogram = [ 0 ] * (buckets_count+1)
    self.incorrect_histogram = [ 0 ] * (buckets_count+1)
    # buckets for levels of bias in non-variations
    self.error_histogram_no_variation = [ 0 ] * (buckets_count+1)

    self.stats = { 'unmapped': 0, 'mapped': 0, 'reference_percent_total': 0., 'error_percent_total': 0., 'incorrect_percent_total': 0., 'total': 0, 'mean_reference': 0, 'mean_error': 0, 'mapped_no_variation': 0, 'unmapped_no_variation': 0 }
    while pos < candidate_fasta.length:
      reference_base = reference.base_at( pos )
      if reference_base is None:
        break # reference is done
      counts = candidate_fasta.count( pos )
      total = sum( [ counts[base] for base in counts ] )

      # TODO any variation
      if pos in donor_vcf.snp_map: # there's a snp here
        self.stats['total'] += 1 # total variations
        if total > 0:
          self.stats['mapped'] += 1
          donor_base = donor_vcf.snp_list[donor_vcf.snp_map[pos]]['alt']
          reference_percent = 1. * counts[reference_base] / total
          error_percent = 1. * ( total - counts[reference_base] - counts[donor_base] ) / total
          self.add_variation_to_stats( reference_percent, error_percent, buckets_count )
          #print "pos", pos, "counts", counts, "ref", reference_base, "vcf", donor_vcf.snp_list[donor_vcf.snp_map[pos]], "ref%", reference_percent, "err%", error_percent
        else: # no coverage
          self.stats['unmapped'] += 1 # is this bias?
      else: # the pos does not contain a variation
        if total > 0: # coverage
          self.stats['mapped_no_variation'] += 1
          error_percent = 1. * ( total - counts[reference_base] ) / total
          self.add_no_variation_to_stats( error_percent, buckets_count )
        else: # no coverage
          self.stats['unmapped_no_variation'] += 1 # is this bias?
        
      #if candidate_base != reference_base and candidate_base != 'N':
      #  target_vcf.snp( pos, reference_base, candidate_base ) # mutation
      pos += 1
    self.reference_histogram = [ 1. * r / self.stats['total'] for r in self.reference_histogram ] # convert to %
    self.error_histogram = [ 1. * r / self.stats['total'] for r in self.error_histogram ] # convert to %
    self.stats['mean_reference'] = 1. * self.stats['reference_percent_total'] / self.stats['mapped']
    self.stats['mean_error'] = 1. * self.stats['error_percent_total'] / self.stats['mapped']

    total_no_variation = self.stats['mapped_no_variation'] + self.stats['unmapped_no_variation']
    self.error_histogram_no_variation = [ 1. * r / total_no_variation for r in self.error_histogram_no_variation ] # convert to %

  def add_variation_to_stats( self, reference, error, buckets_count ):
    # track counts
    self.stats['reference_percent_total'] += reference
    self.stats['error_percent_total'] += error
    self.stats['incorrect_percent_total'] += ( reference + error )
    # find reference bucket
    bucket_size = 1. / buckets_count
    reference_bucket = int( math.floor( reference / bucket_size ) )
    #if reference_bucket >= len(self.reference_histogram):
    #  reference_bucket = len(self.reference_histogram) - 1
    error_bucket = int( math.floor( error / bucket_size ) )
    incorrect_bucket = int( math.floor( ( error + reference ) / bucket_size ) )
    #if error_bucket >= len(self.error_histogram):
    #  error_bucket = len(self.error_histogram) - 1
    self.reference_histogram[reference_bucket] += 1
    self.error_histogram[error_bucket] += 1
    self.incorrect_histogram[incorrect_bucket] += 1

  def add_no_variation_to_stats( self, error, buckets_count ):
    bucket_size = 1. / buckets_count
    error_bucket = int( math.floor( error / bucket_size ) )
    self.error_histogram_no_variation[error_bucket] += 1
