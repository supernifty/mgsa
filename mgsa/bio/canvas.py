import ImageDraw
import os

import fasta

class BamDraw (object):
  def __init__(self, size, target):
    self.size = size
    self.target = target
    self.draw = ImageDraw.Draw(target)
    self.pixels = target.load()

  def add_refs( self, ref_names ):
    self.fasta_names = [ os.path.basename(name) for name in ref_names ]
    self.fasta_lengths = [ fasta.FastaLength( fasta.FastaReader( open( name, 'r' ) ) ).length for name in ref_names ]
    for idx, fasta_info in enumerate( zip(self.fasta_names, self.fasta_lengths) ):
      x = ( idx + 1. ) / ( len(ref_names) + 1 )
      height = 0.8 * fasta_info[1] / max(self.fasta_lengths)
      self.box( (0, 0, 0), x, x + 0.05, 0.1, 0.1 + height )
      self.text( x, 0.05, '%s' % ( fasta_info[0], ) )
      self.text( x, 0.07, '%i' % ( fasta_info[1], ) )

  def text(self, xs, ys, what):
    self.draw.text((int(xs * self.size[0]),int(ys * self.size[1])), what, (0,0,255))

  def box( self, color, xs, xe, ys, ye ):
    xs = int( xs * self.size[0] )
    xe = int( xe * self.size[0] )
    ys = int( ys * self.size[1] )
    ye = int( ye * self.size[1] )
    for x in xrange(xs, xe+1):
      self.pixels[x, ys] = color
      self.pixels[x, ye] = color
    for y in xrange(ys, ye+1):
      self.pixels[xs, y] = color
      self.pixels[xe, y] = color

