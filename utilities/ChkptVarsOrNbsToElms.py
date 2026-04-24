#!/usr/bin/env python

# ChkptVarsOrNbsToElms.py
# Copyright (C) 2026 Wolfgang Tichy
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
from __future__ import print_function

#import numpy as np
import struct
import argparse
import io

# use pythons arg parser
parser = argparse.ArgumentParser(description=
    '''Recreate elms.txt from variables.bin''',
    epilog='''Example:
ChkptVarsOrNbsToElms.py --np 729:9:0 --np 5832:18:1 --OutBoundPats 31,32,33,34,35,36:2 variables.bin''')
parser.add_argument('--np', metavar='POINTS', dest='points', action='append',
        required=True, help='gridpoints in an elm. E.g. --np 729:9:0 prints elms with 729 points as 9*9*9 with pt_typ=0')
parser.add_argument('--OutBoundPats', metavar='OUTBOUNDPATS',
        dest='OutBoundPats',
        help='boundary pats and their pt_typ at outer boundary')
parser.add_argument('file', help='filename')

args = parser.parse_args()


#############################################################################
# read --np args into dict
def points_args_to_dict(points):
  dict = {}
  for np in points:
    sl = np.split(':')
    #print(points, np, sl)
    dict[int(sl[0])] = [int(sl[1]), int(sl[2])]
  return dict

# read --OutBoundPats arg into
def OutBoundPats_arg_to_dict(OutBoundPats):
  dict = {}
  if OutBoundPats != None:
    sl = OutBoundPats.split(':')
    if len(sl) > 0:
      pl = sl[0].split(',')
      pt_typ = int(sl[1])
      for p in pl:
        dict[p] = pt_typ
  return dict


# load data
def load_data(filename, pt_dict, ob_dict):
  buffer = io.StringIO()
  format = 'd'
  size = struct.calcsize(format)
  n_def = 0
  pt_typ_def = 0
  # number of elms will be counted in here
  nelms = 0
  with open(filename, 'rb') as f:
    # read up to 'variable-list-data:\n'
    while True:
      line = f.readline()
      if not line:
        return -1, buffer
      if line == b'variable-list-data:\n':
        break
    # read up to 'variable-list-data:\n'
    while True:
      write_n = 0
      write_pt_typ = 0
      line = f.readline()
      if not line:
        #print('# is this the file end?', file=buffer)
        return nelms, buffer
      if line == b'{\n':
        elmname = f.readline().decode('ascii').strip()
        nelms = nelms + 1 # we found an elm
        p_loc = elmname.split('_')
        pat = p_loc[0] # ob_dict contains str keys
        loc = p_loc[1]
        np = int(f.readline())
        # print what we found
        print(elmname, file=buffer)
        n = pt_dict[np][0]
        pt_typ = pt_dict[np][1]
        # check if pt_typ should be modified on outer boundary
        if pat in ob_dict:
          #print('pat is in ob_dict', file=buffer)
          # check if elm is on outer boundary
          all_odd = True
          for c in loc:
            if int(c) % 2 == 0:
              all_odd = False
              break
          if all_odd:
            #print('cheb', file=buffer)
            pt_typ = ob_dict[pat]
        # check what we print
        if n != n_def:
          write_n = True
        if pt_typ != pt_typ_def:
          write_n = True
          write_pt_typ = True
        # now print
        for i in range(3):
          if write_n:
            print(n, end='', file=buffer)  #fprintf(fp, "%d", n[d]);
          if write_pt_typ:
            print(' ', pt_typ, sep='', end='', file=buffer) #fprintf(fp, " %d", pt_typ[d]);
          if write_n or write_pt_typ:
            print('', file=buffer) # fprintf(fp, "\n");
        # reset defaults
        n_def = n
        pt_typ_def = pt_typ
        print('', file=buffer)
        while True:
          nline = f.readline()
          if not nline:
            #print('# found nothing', file=buffer)
            return -1, buffer
          s = nline.decode('ascii').strip()
          #print('s =', s, file=buffer)
          if s == '}':
            #print('# found }', file=buffer)
            break
          # skip binary var data
          vli_nar = s
          sl = vli_nar.split()
          if len(sl) == 1:
            # skip var data using np
            nskip = np
          elif len(sl) == 4:
            nskip = int(sl[1]) * int(sl[2]) * int(sl[3])
          else:
            #print('# vli_nar is strange:', file=buffer)
            print(vli_nar, file=buffer)
            return -1, buffer
          # skip var data + 1 byte for \n
          #print('# skipping with nskip =', nskip, file=buffer)
          f.seek(nskip*size + 1, 1) #the 1 means seek from current position
  return nelms, buffer

# print
def print_elms(file, pt_dict, ob_dict):
  nelms, buffer = load_data(file, pt_dict, ob_dict)
  # get string from nuffer
  long_string = buffer.getvalue()

  print('number of elms, followed by all elms, their n[3] and optionally their pt_typ[3]')
  print('')
  print('nelms =', nelms)
  print('')
  print(long_string, end='')
  #print(f"{long_string}")
  buffer.close() # Clean up


#############################################################################

# get args
file = args.file
pt_dict = points_args_to_dict(args.points)
ob_dict = OutBoundPats_arg_to_dict(args.OutBoundPats)

#print(pt_dict)
#print(ob_dict)

# load data and then print elms.txt
print_elms(file, pt_dict, ob_dict)
