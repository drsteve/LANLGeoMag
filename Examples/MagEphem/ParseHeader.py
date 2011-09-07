#!/usr/bin/python

import spacepy
import re
import json

filename = '19760724_1976-059.mod.txt'
with open(filename, 'r') as f:
    Lines = f.read()

# isolate header
p = re.compile( r"^#(.*)$", re.M )
h = re.findall( p, Lines )
Header = "".join(h)
#print Header

# isolate JSON field
#s = re.search( r'\s*begin\s+JSON\s*(.*)\s*end\s+JSON', Header )
s = re.search( r'\{\s*(.*)\s*\}', Header )
j = "{" + s.group(1) + "}"
#print j


config_dict = json.loads( j )
spacepy.dictree( config_dict, levels=2, verbose=True )
#print config_dict



