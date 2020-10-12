"""
    Usage: python3 gather_data.py output.json {solver}_data({D:f},{chi:f},{k_d:f}).json ...
"""

import sys
import json
import parse

data = {}

for filename in sys.argv[2:]:
    raw_key = parse.parse("{solver}_data({D:f},{chi:f},{k_d:f}).json", filename.split('/')[-1]).named
    key = (raw_key['solver'], (raw_key['D'], raw_key['chi'], raw_key['k_d']))
    with open(filename, 'r') as f:
        data[str(key)] = json.load(f)

with open(sys.argv[1], 'w') as f:
    json.dump(data, f, separators=(',', ':'))
