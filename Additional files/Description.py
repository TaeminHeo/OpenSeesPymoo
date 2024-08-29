#!/usr/bin/env python
# coding: utf-8

# In[ ]:


### Note 1: section [6] in the 'first py file' should be located in a separate file, 
### Note 2: need additional lines for a fundamental period:
# e.g.
for i in range(0, numEigen):
    lamb = eigenValues[i]
    period = 2 * PI / sqrt(lamb)

#### objective functions and constraints ####

[DispIO, DispLS, DispCP,Area02]=TargetDisp(TotalWeight)  # check 'TargetDisp.m' and 'Bi_Linear2.m'

if DispIO > nStory*Height*0.01:
    disp('목표변위 IO NG')
    DispIO = nStory*Height*0.01*0.98

if DispLS > nStory*Height*0.02:
    disp('목표변위 LS NG')
    DispLS = nStory*Height*0.02*0.98


if DispCP > nStory*Height*0.04:
    disp('목표변위  CP NG')
    DispCP = nStory*Height*0.04*0.98

    
[nodeIO, nodeLS, nodeCP, DriftIO, DriftLS, DriftCP, IndexStep]=ReadOutput(DispIO, DispLS, DispCP, nStory) # check ReadOutput.m
tmpCons = DriftLS(1,1)/AcceptanceCriteria(2,1)
# objective function 1
y(1) = y(1) / 30492000
# objective function 2
y(2) = 1/Area02*10^9
ConsValue = EvalConstraint(nCon, nVar, x, IndexStep, Eccu) # check EvalConstraint.m
cons1 = [tmpCons, ConsValue];

for i in range(0,length(cons1)):
    if cons1(i) > 1.0:
        cons1(i) = cons1(i) - 1.0
    else:
        cons1(i) = 0.0
    
cons =  round(cons1*100)/100

