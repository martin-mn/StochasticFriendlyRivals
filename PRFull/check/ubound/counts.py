import numpy as np
N3=3**16
cls=np.fromfile('cls.bin',dtype=np.uint8); assert len(cls)==N3
idx=np.arange(N3,dtype=np.int64)
D=np.zeros((N3,),dtype=np.int8); x=idx.copy()
dig=np.zeros((16,N3),dtype=np.int8)
for j in range(16): dig[j]=x%3; x//=3
dim=(dig==2).sum(0)
np.save('dim.npy',dim.astype(np.int8))
for c,name in ((1,'C'),(2,'A'),(3,'Z')):
    h=np.bincount(dim[cls==c],minlength=17); print(name,h.tolist(),'total',h.sum())
print('C with p0=1:',((cls==1)&(dig[0]==1)).sum(), 'C total', (cls==1).sum())
print('W dim>=8 class C:', ((cls==1)&(dim>=8)).sum())
# census file
L=open('/Users/martin/Documents/CL1/PartnersRivalsM2f/check/st_eff.txt').read().split('\n')
from collections import Counter
cc=Counter(l.split()[1] for l in L if l.strip()); print('census lines by class',cc)
