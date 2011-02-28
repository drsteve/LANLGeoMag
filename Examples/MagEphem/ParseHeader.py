#!/usr/bin/python

import re
import json

filename = 'JSON.txt'
with open(filename, 'r') as f:
    Lines = f.read()

# isolate header
p = re.compile( r"^#(.*)$", re.M )
m = re.findall( p, Lines )
h = re.findall( p, Lines )
Header = "".join(h)

# isolate JSON field
s = re.search( r'\s*begin\s+JSON\s*(.*)\s*end\s+JSON', Header )
j = s.group(1)

config_dict = json.loads( j )
print config_dict









