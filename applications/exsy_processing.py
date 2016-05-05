# coding: utf-8

from __future__ import print_function
from masai.api import *

# ## read and transform raw Bruker SER files
SYM = True # symmetrize

# epxeriments names
datadir = 'ephedrine'
names = ['160308', '160408']

# create a list of sources
sources = {}

# walk into the experiments directory
for name in names:
    d = os.listdir(os.path.join(datadir,name))
    print("\n --> ", name)

    for expno in d:
        path = os.path.join(datadir, name, str(expno))
        source = Bruker(path)

        # set no output for further processing
        source.silent = True

        # if they are EXSY experiments (read ans processs the SER file)
        if source.par.PARMODE > 0 and source.par.PULPROG == 'exchange.lcs':
            print(expno, end=' ')
            Process(source, transform = 'em 1 --k_shifted 1 0')
            Process(source, transform = 'zf 0')
            Process(source, transform = 'ft')
            Process(source, transform = 'cl 0 12')
            ph = 80
            if name=='160408':
                ph = 87
            Process(source, transform = 'pk %.3f 0 --pivot 9.3'%ph)
            Process(source, axis=0, transform='em 1')
            Process(source, axis=0, transform='zf 2')
            Process(source, axis=0, transform='ft')
            Process(source, axis=0, transform='pk 2.5 0 --select cols --pivot 0')
            Process(source, axis=0, transform = 'cl 0 12')

            # store with an tuple index: mixing time, name
            sources[source.par.D[6], name] = source

# sort the experiment by names and mixing times
k = sorted(sources.keys())
k3=[]
k4=[]
for i in k:
    if i[1].endswith('308'):
        k3.append(i)
    else:
        k4.append(i)
k3 = sorted(k3)
k4 = sorted(k4)

# Show control plot
id = -12
s = sources[k3[id]]

row = s.get_row(9.4)
p = Plot1d(row, figsize = (5, 4), ylim=(0,150000), display=False)
col = s.get_col(9.4)
Plot1d(col, axe = p.axe, display=False)

s = sources[k4[id]]

row = s.get_row(9.4)
Plot1d(row, axe = p.axe, display=False)
col = s.get_col(9.4)
Plot1d(col, axe = p.axe, display=False)

source1 = sources[k3[id]]
Plot2d(source1, ylim=(0,12), figsize=(4.5,4.5), proj='sky', display=False)
Plot2d(source1.sym, ylim=(0,12), figsize=(4.5,4.5), proj='sky', display=False)

source2 = sources[k4[id]]
Plot2d(source2, ylim=(0,12), figsize=(4.5,4.5), proj='sky', display=False)
Plot2d(source2.sym, ylim=(0,12), figsize=(4.5,4.5), proj='sky', display=False)

if SYM:
    for k, source in sources.items():
        sources[k] = source.sym


id = -12
s = sources[k3[id]]

row = s.get_projection()
p = Plot1d(row, figsize = (5, 4), display=False)

s = sources[k4[id]]
row = s.get_projection()
Plot1d(row, axe = p.axe, display=False)#
thdelta  = np.array([7.8, 8.0, 7.9, 2.1, 5.4, 0.9, 3.0, 5.3, 4.6, 1.8])+1.8
thdeltadA  = np.array([6.8, 7.3, 7.5, 1.9, 5.3, 1.8, 3.4, 4.5, 8.3, 3.2])+1.8
thdeltadB  = np.array([6.2, 7.4, 7.7, 1.5, 5.4, 1.5, 4.9, 5.5, 3.9, 3.7])+1.8


for i, d in enumerate(thdelta):
    p.axe.axvline(d, ymax=.7, ls=':')
    p.axe.text(d+.1, 300000-i*2000, 'H%d'%(i+1), fontsize=7, color="blue")

for i, d in enumerate(thdeltadA):
    p.axe.axvline(d, ymin= .3, ls=':', color='green')
    p.axe.text(d+.1, 250000-i*2000, 'H%d'%(i+1), fontsize=7, color="green")

for i, d in enumerate(thdeltadB):
    p.axe.axvline(d, ymin= .3, ls=':', color='red')
    p.axe.text(d+.1, 200000-i*2000, 'H%d'%(i+1), fontsize=7, color="red")


# In[ ]:

# lets try to use PCA to find the components
p = Pca(sources)
print(p.quality(), p.npc)
p.model()
# performs an automatic baseline correction of the PC's
p.basecorr()
newsources = p.transform()


pl.figure()
l = p.loadings
pl.plot(l.columns, l.T.values[:,:6])

# In[ ]:

# Attempt with NMF

# we assume homogeneous dataframe same index and columns
d = []
for k, v in newsources.iteritems():
    ir = v.data.values.real
    m=ir.min()
    d.append(ir-m)

x = v.data.columns
X = np.concatenate(d, axis=0)

from sklearn.decomposition import NMF
nc=10
model = NMF(n_components=nc, init='random', random_state=0)
model.fit(X)

PC = model.components_



id = -12
s = newsources[k3[id]]

row = s.get_row(9.4)
p = Plot1d(row, figsize = (5.5, 4), display=False)

s = newsources[k4[id]]
row = s.get_row(9.4)
Plot1d(row, axe = p.axe, display=False)
p.axe.plot(x, PC[:nc].T * 50, lw=.3)


def adjustzones(newsources, zones):
    idzones= []
    for iz, zone in enumerate(zones):
        index = newsources[newsources.keys()[0]].data.columns.values
        select = np.argwhere((index>zone[0]) & (index<zone[-1])).flatten()
        i1, i2 = select[0], select[-1]
        pos1 = index[i1]
        pos2 = index[i2]
        if i2 < i1:
            i1, i2 = i2, i1
            pos1, pos2 = pos2, pos1
        idzones.append((i1,i2))
        zones[iz] = (pos1, pos2)
    return idzones, zones

_zones = [(2.5, 3.0),
          (3.7, 4.1),
          (4.4, 4.7),
          (4.8, 5.2),
          (5.7, 6.1),
          (6.6, 7.1),
          (8.3, 8.8),
          (9.0, 9.3),
          (9.3, 9.6)]
izones, zones = adjustzones(newsources, _zones)

# graphical check
s = newsources[k3[id]]
p = Plot2d(s, proj='sky', figsize=(4.5,4.5), display=False)

for z1 in zones:
    p.add_xproj_zone(z1)
    p.add_yproj_zone(z1)
    for z2 in zones:
        p.add_zone(z1,z2)

s = newsources[k4[id]]
p = Plot2d(s, proj='sky', figsize=(4.5,4.5), display=False)

for z1 in zones:
    p.add_xproj_zone(z1)
    p.add_yproj_zone(z1)
    for z2 in zones:
        p.add_zone(z1,z2)

pl.show()
sys.exit()

# In[ ]:

d = []
d2 = []
y = []
y2 = []

for k in sorted(k3):
    item = newsources[k]
    ir = item.data.values.real
    m=ir.min()
    d.append(ir-m)
    y.append(k[0])

for k in sorted(k4):
    item = newsources[k]
    ir = item.data.values.real
    m=ir.min()
    d2.append(ir-m)
    y2.append(k[0])

c2 = np.dstack(d2)
c2 = c2.swapaxes(0,2)

c = np.dstack(d)
c = c.swapaxes(0,2)

x = item.data.columns

x1 = 0
for inc in range(0,306,18):
    fig = pl.figure()
    x2 = inc
    print(x[x1],x[x2])
    pl.plot(y, (c[:,x1:x2,x1:x2].sum(axis=1)).sum(axis=1), 'o-', ms=5)
    pl.plot(y2, (c2[:,x1:x2,x1:x2].sum(axis=1)).sum(axis=1), 'o-', ms=5)
    x1 = x2

    pl.xlim(0,0.2)
    pl.show()


# In[ ]:

d = []
d2 = []
y = []
y2 = []

for k in sorted(k3):
    item = sources[k]
    ir = item.data.values.real
    m=ir.min()
    d.append(ir-m)
    y.append(k[0])

for k in sorted(k4):
    item = sources[k]
    ir = item.data.values.real
    m=ir.min()
    d2.append(ir-m)
    y2.append(k[0])

c2 = np.dstack(d2)
c2 = c2.swapaxes(0,2)

c = np.dstack(d)
c = c.swapaxes(0,2)

x = item.data.columns

x1 = 0
for inc in range(0,306,18):
    fig = pl.figure()
    x2 = inc
    print(x[x1],x[x2])
    pl.plot(y, (c[:,x1:x2,x1:x2].sum(axis=1)).sum(axis=1), 'o-', ms=5)
    pl.plot(y2, (c2[:,x1:x2,x1:x2].sum(axis=1)).sum(axis=1), 'o-', ms=5)
    x1 = x2

    pl.xlim(0,0.2)
    pl.show()


# In[ ]:


