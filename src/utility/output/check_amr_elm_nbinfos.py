#!/usr/bin/env python

import argparse

# use pythons arg parser
parser = argparse.ArgumentParser(description=
    '''Check if amr_elm_nbinfo contains the correct eid's.
       This assumes the nmesh-function write_amr_elm_nbinfos has been
       called to create the amr_elm_nbinfos.txt file.''',
    epilog='''Example:
    check_amr_elm_nbinfos MPA1/amr_elm_nbinfos.txt''')
parser.add_argument('file', help='pathname of amr_elm_nbinfos file')

args = parser.parse_args()

# open file bak_list.txt
#print("# Opening:", args.file)
with open(args.file, 'r') as f:
    lines = f.readlines()

# dict block contains the blocks with elms from file
block = {}
for line in lines:
    #print('line =', line)
    # skip empty lines
    if len(line) == 0:
        continue
    # if 1st char is space we have a nb line that we will skip
    if line[0] == ' ':
        continue
    # if 1st char is not a number we have a header
    if (not line[0].isdigit()):
        header = line
        elms = {} # start new dict
        block[header] = elms
        continue
    # if 1st char is a number we have an elm
    if line[0].isdigit():
        #print(line[0], line[0].isdigit())
        elm = line.split(':')[0]
        eid = line.split(' ')[1]
        elms = block[header]
        elms[elm] = eid

#print(block)

# go through lines again and check every nb eid
for line in lines:
    # skip empty lines
    if len(line) == 0:
        continue
    # if 1st char neither a number nor space we have a header
    if not (line[0].isdigit() or line[0] == ' '):
        header = line
        continue
    # if 1st char is a number we have an elm
    if line[0].isdigit():
        elm = line.split(':')[0]
    # if 1st char is space we have a nb line that we will check now
    if line[0] == ' ':
        #print(line)
        face = line.split(':')[0].strip()
        pos = line.find(':') + 1
        pos = line.find(':', pos) + 1
        nbstr = line[pos:].strip()
        #print(pos, nbstr)
        nbs = nbstr.split(' ')
        for nb in nbs:
            nbelm = nb.split(':')[0]
            nbeid = nb.split(':')[1]
            elms = block[header]
            #print(nbelm, nbeid)
            #print(elms)
            #print('nb', nb, 'in', elm, elms[nbelm])
            if elms[nbelm] != nbeid:
                print('block',header.strip(), 'elm',elm, ': nb',face,nb,
                      'HAS WRONG eid!!!')
