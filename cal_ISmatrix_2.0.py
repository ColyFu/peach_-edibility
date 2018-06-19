#!/usr/bin/env python
import sys
import os
import re
import multiprocessing
from collections import defaultdict
import random
import time

Start = time.time()

iscore = defaultdict(lambda: 0)
iscorec = defaultdict(lambda: 0)
iscorev = defaultdict(lambda: 0)
samples = []
q = multiprocessing.Queue()
iscorev = {}
cpu = 8


for line in open(sys.argv[2]):
    if line.startswith('#'):
        ff = open(sys.argv[1], 'r')
        ff.seek(int(line.strip('\n').split('\t')[1]), 0)
        samples = ff.readline().strip('\n').split('\t')[9:]
    else:
        q.put(line.strip('\n'), block=False)


def calindex(line, samples):
    index = line.split('\t')
    ff = open(sys.argv[1], 'r')
    ff.seek(int(index[1]), 0)
    iscore_t = defaultdict(lambda: 0)
    iscorec_t = defaultdict(lambda: 0)
    samnum = len(samples)
    for line in ff:
        info = line.strip('\n').split('\t')
        if info[0] != index[0]:
            break
        fmt = info[8].split(':').index('AD')
        fall = []
        for each in info[9:]:
            try:
                ref, alt = each.split(':')[fmt].split(',')
                fall.append(float(ref) / (int(ref) + int(alt)))
            except:
                fall.append('NA')

        for i in xrange(samnum):
            for j in xrange(i + 1, samnum):
                key = samples[i] + ':' + samples[j]
                try:
                    iscore_t[key] += 1 - abs(fall[i] - fall[j])
                    iscorec_t[key] += 1
                except:
                    continue
    rdm = ''.join([str(int(random.random() * 10)) for i in xrange(8)])
    ff = open('IS.cal.temp' + rdm, 'w')
    for each in iscore_t:
        ff.write("%s\t%s\t%s\n" % (each, iscore_t[each], iscorec_t[each]))


def consumer(q, samples):
    while True:
        try:
            line = q.get(block=False)
            calindex(line, samples)
        except:
            break


con = {}
for each in xrange(cpu):
    con[each] = multiprocessing.Process(target=consumer, args=(q, samples))
    con[each].daemon = True
    con[each].start()

for each in con:
    con[each].join()

ss = len(samples)
tmpf = os.listdir('./')
for each in tmpf:
    if re.search('IS.cal.temp\d{6}', each):
        for line in open(each):
            info = line.strip('\n').split('\t')
            iscore[info[0]] += float(info[1])
            iscorec[info[0]] += float(info[2])

os.system("/bin/rm IS.cal.temp*")

for i in xrange(ss):
    for j in xrange(i + 1, ss):
        key = samples[i] + ':' + samples[j]
        iscorev[key] = iscore[key] / iscorec[key]
        iscorev[samples[j] + ':' + samples[i]] = iscorev[key]

fout = open(sys.argv[3], 'w')
fout.write("Sample" + '\t' + '\t'.join(samples) + '\n')
for i in samples:
    line = i
    for j in samples:
        if i == j:
            line += '\t1'
        else:
            line += '\t%s' % (iscorev[i + ':' + j])
    fout.write(line + '\n')

End = time.time()

Cons = End - Start
m, s = divmod(Cons, 60)
h, m = divmod(m, 60)

print "%d:%02d:%02d be used" % (h, m, s)
