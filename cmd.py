
import argparse
import mgsa

class CommandLine(object):
  def update(self, complete, what):
    print complete, what
    return True

  def align(self, genomes, sequences):
    aligner = mgsa.MultipleGenomeSequenceAligner(args.genomes)
    aligner.align( args.sequences, self )

parser = argparse.ArgumentParser(description='Assemble sequence reads against multiple genomes')
parser.add_argument('-s', '--sequences', nargs=2, type=file, required=True,
                   help='Paired end FASTQ sequence files')
parser.add_argument('-g', '--genomes', nargs='+', type=file, required=True,
                   help='Reference genome FASTA file(s)')

args = parser.parse_args()

processor = CommandLine()
processor.align( args.genomes, args.sequences )
