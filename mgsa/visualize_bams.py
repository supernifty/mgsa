
import argparse
import Image
import sys

import bio
import config

IMAGE_SIZE = (1200, 800)

parser = argparse.ArgumentParser(description='Visualize BAMs')
parser.add_argument('--bams', dest='bams', nargs='+', help='bam files to analyze')
parser.add_argument('--refs', dest='refs', nargs='+', help='fasta files to analyze')
args = parser.parse_args()

draw = bio.BamDraw( IMAGE_SIZE )
draw.add_refs( args.refs )
draw.target.show()

