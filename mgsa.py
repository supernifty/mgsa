import threading
import datetime
import time

class MultipleGenomeSequenceAlignerThread(threading.Thread):
  '''
    Threaded wrapper to MultipleGenomeSequenceAligner
  '''

  def __init__(self, genomes, sequence):
    threading.Thread.__init__(self)
    self.worker = MultipleGenomeSequenceAligner(genomes)
    self.sequence = sequence
    self.status = []
    self.running = True
    self.start()

  def run(self):
    self.worker.align(self.sequence, self)
    self.running = False

  def stop_align(self):
    self.running = False

  def update( self, completion, status ):
    self.status.append( { 'when': datetime.datetime.now().strftime("%Y%m%d-%H%M%S"), 'what': status, 'complete': completion } )
    return self.running

  def clear_status( self ):
    result = self.status[:]
    self.status = []
    return result

class InterruptedError(Exception):
  pass

class MultipleGenomeSequenceAligner(object):
  '''
    Main class that accepts multiple genomes and attempts to align a sequence to them
  '''

  def __init__(self, genomes):
    self.genomes = genomes
    self.status = None

  def _update_status(self, completed, what):
    if self.status is not None:
      result = self.status.update( completed, what )
      if not result:
        raise InterruptedError

  def align(self, sequence, status):
    self.status = status
    try:
      self._update_status( 0, "starting..." ) 
      time.sleep(5)
      self._update_status( 100, "finished" )
    except InterruptedError:
      pass
