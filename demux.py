#!/usr/bin/env python
# -*- coding:utf-8 -*-

import seeq
import sys
import gzip

from gzopen import gzopen

IND = [dict(),dict()]
params = {}
lines = ['']*4

# Helper functions.
def getindex(matcher, txt):
   suffix = matcher.matchSuffix(txt, False) or ''
   return suffix[:4]

def main(fname):
   lineno = 0
   mode = 2
   
   # Parse parameter file
   cf = open(fname)
   for line in cf:
      lineno += 1
      line = ''.join(line.split())
      if line == '' or line[0] == '#':
         continue
      elif len(line.split('=')) == 2:
         [param,value] = line.split('=')
         params[param] = value
      elif len(line.split(':')) == 2:
         section = line.split(':')[0]
         if section == 'dna-index':
            mode = 0
         elif section == 'rna-index':
            mode = 1
         else:
            print "error in parameter file: Uknown section '{}' in {}, line {}.".format(section,fname,lineno)
            sys.exit(1)
      elif mode < 2 and len(line.split(',')) == 2:
         [index,fout] = line.split(',')
         if IND[mode].has_key(index):
            print "duplicate index '{}' in {}, line {}".format(index,fname,lineno)
            sys.exit(1)
         IND[mode][index] = fout
      else:
         print "unknown parameter '{}' in {}, line {}".format(line,fname,lineno)
         sys.exit(1)

   # Check parameters
   error = 0
   if not params.has_key('bfs'):
      print "missing parameter in {}: 'bfs' must be defined.".format(fname)
      error = 1
   if not params.has_key('dist'):
      print "missing parameter in {}: 'dist' must be defined.".format(fname)
      error = 1
   if not params.has_key('dna-seqfile'):
      print "missing parameter in {}: 'dna-seqfile' must be defined.".format(fname)
      error = 1
   if not params.has_key('rna-seqfile'):
      print "missing parameter in {}: 'rna-seqfile' must be defined.".format(fname)
      error = 1
   if error:
      sys.exit(1)

   # Compile flanking sequence
   T7 = seeq.compile(params['bfs'], int(params['dist']))

   FDICT = [dict(),dict()]
   for index,fname in IND[0].items():
      if not FDICT[0].has_key(fname):
         FDICT[0][fname] = gzip.open(fname, 'wb')
         
   for index,fname in IND[1].items():
      if not FDICT[1].has_key(fname):
         FDICT[1][fname] = gzip.open(fname, 'wb')

   try:
      # Demultiplex DNA indices
      with gzopen(params['dna-seqfile']) as f:
         # Read fastq file.
         for lineno,line in enumerate(f):
            lines[lineno%4] = line
            if lineno % 4 == 3:
               index = getindex(T7, lines[1])
               if IND[0].has_key(index):
                  f = FDICT[0][IND[0][index]]
                  for l in lines:
                     f.write(l)

      # Demultiplex RNA indices
      with gzopen(params['rna-seqfile']) as f:
         # Read fastq file.
         for lineno,line in enumerate(f):
            lines[lineno%4] = line
            if lineno % 4 == 3:
               index = getindex(T7, lines[1])
               if IND[1].has_key(index):
                  f = FDICT[1][IND[1][index]]
                  for l in lines:
                     f.write(l)
   finally:
      for f in FDICT[0].values(): f.close()
      for f in FDICT[1].values(): f.close()
      
                  
if __name__ == '__main__':
   main(sys.argv[1])

