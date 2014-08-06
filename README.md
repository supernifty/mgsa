
Multiple Genome Sequence Aligner
================================
This work is in progress and not ready for public consumption.
Trying to use this code is strongly discouraged.

Setup
-----
* config.py - edit this to paths appropriate for your system

Included Files
--------------
* cmd.py - command line based interface
* mgsa.py - multi genome sequence aligner
* bio.py - helper classes
* ui.py - web based user interface
* static/ - included static files for the web based UI
* state/ - details of current web based session

Example Usage
-------------
python cmd.py -s ../data/tiny_r1.fastq ../data/tiny_r2.fastq -g ../data/tiny_reference.fasta

