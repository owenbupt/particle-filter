#!/usr/bin/env python
# coding=utf-8
# Author = huhuhushan

import math
import random
x = 0.1
Q = 1
R = 1
tf = 50
N = 100

xhat = x
P = 2
xhatPart = x
xpartminus = []

xpart = [x + math.sqrt(P) * random.normalvariate(0, 1) for i in range(N)] # 
# print(xpart)
q = [] 
xArr = [x]
xhatPartArr = [xhatPart]
for k in range(tf):
    xpartminus = []
    q = [] 
    x = 0.5 * x + 25 * x / (1 + pow(x, 2)) + 8 * math.cos(1.2*(k-1)) + math.sqrt(Q) *random.normalvariate(0, 1) 

    y = math.pow(x, 2) / 20 + math.sqrt(R) * random.normalvariate(0, 1) 
    for i in range(N):
        xpartminus.append(0.5 * xpart[i] + 25 * xpart[i] / (1 + math.pow(xpart[i], 2)) + 8 * math.cos(1.2 * (k - 1)) + math.sqrt(Q) * random.normalvariate(0, 1))

        ypart = math.pow(xpartminus[i], 2) / 20
        vhat = y - ypart
        q.append(1 / math.sqrt(R) / math.sqrt(2 * math.pi) * math.exp(-math.pow(vhat, 2) /2 / R))


    qsum = 0
    

    for i in range(N):
        qsum += q[i]

    q = [q[i] / qsum for i in range(N)]

    for i in range(N):
        u = random.random()
        qtemsum = 0
        for j in range(N):
            qtemsum += q[j]
            if qtemsum >= u:
                xpart[i] = xpartminus[j]
                break
    xpart_sum = 0
    for i in range(N):
        xpart_sum += xpart[i]
    xhatPart = xpart_sum / N
    xArr.append(x)
    xhatPartArr.append(xhatPart)

print([xhatPartArr[i] - xArr[i] for i in range(tf)])


