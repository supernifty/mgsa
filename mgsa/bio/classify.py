
# train and classify variations

class ReadDepthLogisticClassifier (object):
  '''
    build a classifier using local read depth data with 2 possible labels
  '''

  def __init__( self, fasta, read_length ):
    '''
      @fasta: ProbabilisticFasta     
      @read_length: determines how local the features will be (i.e. +/- read_legnth)
    '''
    self.fasta = fasta
    self.read_length = read_length
    self.window_start = -read_length
  
  def add_instance( pos, positive ):
    '''
      @pos: position on fasta
      @positive: is this a positive example
    '''
    pass

  def classify( pos ):
    '''
      return a probability of positive classification
    '''
    pass
