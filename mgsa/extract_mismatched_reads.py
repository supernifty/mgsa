
# given input from compare_bams, generate a sam file with matching mismatches based on specified criteria
# e.g. python extract_mismatched_reads.py --min_distance 50 --mappable True ../../data/e-coli-mg1655_snp_ins_del_50_x1_bwa.bam < ../../data/compare_bams3.log > mappable_test.sam

import argparse
import sys

import bio
import config

parser = argparse.ArgumentParser(description='Find mismatched reads matching criteria')
parser.add_argument('bam', help='bam files to extract reads')
parser.add_argument('--min_distance', dest='min_distance', type=int, default=0, help='only consider alignments with at least this distance from the true location')
parser.add_argument('--max_distance', dest='max_distance', type=int, default=1e12, help='only consider alignments with at most this distance from the true location')
parser.add_argument('--mappable', dest='mappable', type=bool, default=False, help='only consider alignments that have a mappable location')

args = parser.parse_args()

allowed_tags = set()
for line in sys.stdin:
  fields = line.strip().split(',')
  if len(fields) == 5 and fields[0].isdigit():
    distance = abs(int(fields[3]))
    mappable = int(fields[4]) == 0
    if distance >= args.min_distance and distance <= args.max_distance:
      if mappable and args.mappable or not mappable and not args.mappable:
        allowed_tags.add( fields[2] )
bio.log_stderr( '%i allowed tags' % len(allowed_tags) )

bio.SamFilter( sam_fh=bio.BamReaderExternal( config.BAM_TO_SAM, args.bam ), target_fh=sys.stdout, allowed_tags=allowed_tags, log=bio.log_stderr )
