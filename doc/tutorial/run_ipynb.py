# Here we use runipy to genrate all output from the tutorial notebooks
# and may be detect  exceptions
import os
from runipy.notebook_runner import NotebookRunner
from IPython.nbformat.current import read, write

ldir = os.listdir('.')
print(ldir)
for item in ldir:
    if os.path.splitext(item)[1]==".ipynb":
        ipynb = os.path.join(item)
        print "\n >>>>>> running %s" % ipynb
    else:
        continue
    notebook = read(open(ipynb), 'json')
    r = NotebookRunner(notebook)
    print(r.run_notebook())
    print('******')
    write(r.nb, open(ipynb, 'w'), 'json')