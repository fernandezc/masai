"""This listing is gonna be uploaded to the server
"""
import os
import sys
args = sys.argv
path = args[1]
catalog = {}
count = 0
for root, dirs, files in os.walk(path):
    for fic in files:
        if fic in ['ser','fid']:
            name, expno = os.path.split(root)
            if name not in catalog.keys():
                catalog[name]=[]
            # add some information
            catalog[name].append(expno)
            count += 1
            
if count == 0:
    print 'no file found'

print catalog