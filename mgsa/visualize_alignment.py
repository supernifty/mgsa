
import cgi
import datetime
import BaseHTTPServer
import SocketServer
import webbrowser
import json
import re

import bio

PORT = 8421
STATE = 'state/state.json'

class Session(object):
  processor = None # class variable

  def __init__(self):
    # get current state
    try:
      f = open(STATE, 'r')
      self.state = json.loads( f.read() )
      f.close()
    except IOError:
      self.state = {}
 
  def save(self):
    f = open(STATE, 'w')
    f.write( json.dumps(self.state) )
    f.close()

  def log(self, msg):
    if 'log' not in self.state:
      self.state['log'] = []
    when = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    self.state['log'].append( { 'when': when, 'what': msg } )

  def add_genome(self, fh, title):
    if 'genomes' not in self.state:
      self.state['genomes'] = []
    # save data
    location = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    f = open( 'state/%s.ref' % location, 'w' )
    f.write( fh.read() ) 
    # remember
    self.state['genomes'].append( {'title': title, 'location': location } )
    # log
    self.log( 'added genome "%s" as "%s"' % ( title, location ) )

  def add_alignment(self, fh, title):
    if 'alignments' not in self.state:
      self.state['alignments'] = []
    # save data
    location = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    f = open( 'state/%s.sam' % location, 'w' )
    f.write( fh.read() ) 
    # remember
    self.state['alignments'].append( {'title': title, 'location': location } )
    # log
    self.log( 'added alignment "%s" as "%s"' % ( title, location ) )

  def update_status(self):
    self.save()

class Visualizer(object):
  def __init__(self, ref, sam, start, end):
    self.ref_fh = open(ref)
    self.sam_fh = open(sam)
    reader = bio.FastaReader(self.ref_fh)
    current = 0
    self.ref_data = ''
    for line in reader.items():
      if current >= start and current <= end: # take it all
        self.ref_data += line
      elif current < start and current <= end: # take partial
        self.ref_data += line[start-current:]
      current += len(line)
      if current > end:
        break
    if len(self.ref_data) > end - start + 1:
      self.ref_data = self.ref_data[:end - start + 1]
    
  def visualize(self):
    return { 'ref': self.ref_data, 'sam': '' }
    
class AlignmentVisualizerHTTPRequestHandler(BaseHTTPServer.BaseHTTPRequestHandler):
  def do_HEAD(self):
    self.send_response(200)
    self.send_header("Content-type", "text/html")
    self.end_headers()

  def do_GET(self):
    """Respond to a GET request."""
    self.send_response(200)
    while self.path.startswith( '/' ):
      self.path = self.path[1:]
    if self.path == '':
      self.path = 'visualize.htm'
    if self.path.find( '..' ) == -1:
      f = open( 'static/%s' % self.path, 'r' )
      if self.path.endswith( '.htm' ):
        self.send_header("Content-type", "text/html")
      if self.path.endswith( '.css' ):
        self.send_header("Content-type", "text/css")
      if self.path.endswith( '.js' ):
        self.send_header("Content-type", "application/js")
      self.end_headers()
      self.wfile.write( f.read() )

  def do_POST(self):
    """Respond to a POST request."""
    session = Session()
   
    while self.path.startswith( '/' ):
      self.path = self.path[1:]
    
    if self.path == 'status':
      self.send_response(200)
      self.send_header("Content-type", "application/json")
      self.end_headers()
      session.update_status()
      self.wfile.write( json.dumps( session.state ) )
      return
      
    if self.path == 'add_genome':
      form = cgi.FieldStorage(fp=self.rfile, headers=self.headers, environ={'REQUEST_METHOD':'POST', 'CONTENT_TYPE':self.headers['Content-Type'], })
      upload = form['upload']

      # extract basename of input filename, remove non-alphanumeric characters
      title = re.split( '[\\\/]', upload.filename )[-1]
      session.add_genome( upload.file, title )
      
    if self.path == 'add_alignment':
      form = cgi.FieldStorage(fp=self.rfile, headers=self.headers, environ={'REQUEST_METHOD':'POST', 'CONTENT_TYPE':self.headers['Content-Type'], })
      upload = form['upload']

      # extract basename of input filename, remove non-alphanumeric characters
      title = re.split( '[\\\/]', upload.filename )[-1]
      session.add_alignment( upload.file, title )
      
    if self.path == 'visualize':
      self.send_response(200)
      self.send_header("Content-type", "application/json")
      self.end_headers()
      length = int(self.headers.getheader('content-length'))
      params = cgi.parse_qs(self.rfile.read(length), keep_blank_values=1)
      result = Visualizer( 'state/%s.ref' % params['ref'][0], 'state/%s.sam' % params['sam'][0], int(params['start'][0]), int(params['end'][0]) ).visualize()
      self.wfile.write( json.dumps( result ) )
      return

    # save state
    session.save()

    self.send_response(302)
    self.send_header("Location", "/")
    self.end_headers()

handler = AlignmentVisualizerHTTPRequestHandler
httpd = SocketServer.TCPServer(("", PORT), handler)

print "serving at port", PORT
webbrowser.open( "http://127.0.0.1:%i/" % PORT )
httpd.serve_forever()

