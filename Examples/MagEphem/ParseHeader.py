#!/usr/bin/python

import spacepy
import re
import json

filename = 'JSON.txt'
with open(filename, 'r') as f:
    Lines = f.read()

# isolate header
p = re.compile( r"^#(.*)$", re.M )
#m = re.findall( p, Lines )
h = re.findall( p, Lines )
Header = "".join(h)
print Header
exit()

# isolate JSON field
#s = re.search( r'\s*begin\s+JSON\s*(.*)\s*end\s+JSON', Header )
s = re.search( r'\{(.*)\}', Header )
j = s.group(1)


config_dict = json.loads( j )
spacepy.dictree( config_dict, levels=2, verbose=True )
#print config_dict



