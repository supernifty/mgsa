
Multiple Genome Sequence Aligner
================================
[![Build Status](https://travis-ci.org/supernifty/mgsa.svg?branch=master)](https://travis-ci.org/supernifty/mgsa)
[![Coverage Status](https://coveralls.io/repos/supernifty/mgsa/badge.png?branch=master)](https://coveralls.io/r/supernifty/mgsa?branch=master)

This work is in progress and not ready for public consumption.
Trying to use this code is strongly discouraged.

Setup
-----
* config.py - edit this to paths appropriate for your system

Dependencies
------------

Tests
-----
Run test.sh

Included Files
--------------
* __init__.py      
* batch_plot.py - some example plotting routines          
* build_repeated_fasta.py
* evaluate.py
* evaluate_reads.py
* find_repeat_probability.py
* find_repeats.py
* find_repeats_large.py
* generate_all_reports.py - builds all plots used in the final report
* generate_consensus.py  
* generate_mutation.py   
* generate_reads.py      
* generate_vcf.py  
* helpers.py       
* main.py             
* mapper.py - a naive, simple aligner
* mapper_selector.py - runs an aligner
* mgsa.py - multi genome sequence aligner
* pipeline_batch.py - start a batch of tests specified by the provided batch config file
* pipeline_test.py   
* sam_to_vcf.py
* ui.py - web based user interface
* wgsim_pipeline.py

* batch/ - batch files that run multiple alignment tests
* bio/ - common bio classes
  * bio.py
  * fasta.py
  * __init__.py
  * report.py
  * sam.py
  * vcf.py

* out/ - results of alignment tests
* state/ - details of current web based session
* static/ - included static files for the web based UI
* tests/

Deprecated
----------
* cmd.py - command line based interface

Pipeline Config Files
---------------------

Example Usage
-------------

Pipeline
* python pipeline_batch.py batch/pipeline_batch_bias_ecoli_2.cfg out/pipeline_batch_bias_ecoli_2.cfg 

Deprecated
----------

python cmd.py -s ../data/tiny_r1.fastq ../data/tiny_r2.fastq -g ../data/tiny_reference.fasta

