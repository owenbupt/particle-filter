#!/usr/bin/env python
# coding=utf-8
# Author = huhuhushan

import math
import random
x1 = 2500
y1 = 25000
vx1 = 0
vy1 = -222
ax1 = 7
ay1 = 7
Ts = 2
tf = 50
N = 100

xhat = x
P = 2
xhatPart = x
xpartminus = []

x1part = [x1 for i in range(N)]
y1part = [y1 for i in range(N)]
vx1part = [vx1 for i in range(N)]
vy1part = [vy1 for i in range(N)]
# print(xpart)
q = [] 
xArr = [x]
xhatPartArr = [xhatPart]
for k in range(tf):
    x1partminus = []
    y1partminus = []
    vx1partminus = []
    vy1partminus = []
    q = [] 
    x1 = x1 + Ts * vx1 + Ts ** 2 / 2 * ax1 * random.normalvariate(0, 1)
    y1 = y1 + Ts * vy1 + Ts ** 2 / 2 * ay1 * random.normalvariate(0, 1)
    vx1 = vx1 + Ts * ax1 * random.normalvariate(0, 1)
    vy1 = vy1 + Ts * ay1 * random.normalvariate(0, 1)
    for i in range(N):
        x1partminus.append(x1part[i] + )
        x1 = x1 + Ts * vx1 + Ts ** 2 / 2 * ax1 * random.normalvariate(0, 1)
        y1 = y1 + Ts * vy1 + Ts ** 2 / 2 * ay1 * random.normalvariate(0, 1)
        vx1 = vx1 + Ts * ax1 * random.normalvariate(0, 1)
        vy1 = vy1 + Ts * ay1 * random.normalvariate(0, 1)
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


