## Install 
```shell
sh autogen.sh
./configure
make -j 10 && make install

pip install cython --upgrade --user # cython >= 0.28.2 is ok
python setup.py build_ext --inplace
```

## Run
```py
from sv_helper import LocalSGAWrapper
sga=LocalSGAWrapper('test', 0.01, 30, 75)
sga.run(['GGCAGACAGATCACTTGAGGTCAGAAATTGAAGGCCAGCCTGGCCAACATGGTGAAACCCCATCTCTACTAAAAA', 'ATTGAAGGCCAGCCTGGCCAACATGGTGAAACCCCATCTCTACTAAAAATACAAAAATTAGCCAGGCATGCTGGC'],2)

```

### core dump

```py
from sv_helper import LocalSGAWrapper
sga=LocalSGAWrapper('test', 0.01, 30, 75)

from random import randint, random
alphabet_set = 'ATCG'
src = ''.join(alphabet_set[randint(0, 3)] for _ in range(1000))
s = []
for _ in range(100):
	if random() < 0.5:
		st = randint(100,200)
	else:
		st = randint(600,700)
	len=randint(60,80)
	s.append(src[st:st+len])

sga.run(s, 3) # >2 will core dump, 可能和聚类数有关？
```