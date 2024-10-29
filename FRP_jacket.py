import numpy as np
import pandas as pd
import math
from openseespy.opensees import *

def opensees_configure(x, nVar=6, nObj=3, nCon=12+1):

    # Clear model
    wipe()

    # Create ModelBuilder (with two-dimensions and 3 DOF/node)
    model('basic', '-ndm', 2, '-ndf', 3)

    # Set parameters for overall model geometry
    CRB = 0.5
    CRC = 0.7

    nStory = 3
    nSpan = 3

    Height = 3000
    Span = 6000

    EffectiveDepth = 6000/1.0

    # Loads calc
    DL = 5.0*1000/1000000
    LL = 2.4*1000/1000000

    tmp = 1.0*(DL+0.25*LL)
    beamload = (-1)*(DL+0.25*LL)*(1000)
    unitWeight = np.ones((3,1))*(DL)
    ShearPoint = np.zeros((6,4))
    nColumnMaterial = (nStory)*(nSpan+1)
    nColumn = (nStory)*(nSpan+1)
    Ast = np.zeros((nColumnMaterial,1))
    Eccu = np.zeros(((nStory)*(nSpan+1),1))

    g = 9806.65

    # Pushover control
    IDctrlNode = 41
    IDctrlDOF = 1
    nStep = 500
    Dmax = Height*nStory*0.04*1.0
    Dincr = Dmax/nStep

    # Variables for objective functions and constraints
    y = np.zeros((1,nObj))
    cons = np.zeros((1,nCon))

    tf = 0.33
    AcceptanceCriteria = [0.01, 0.02, 0.04]
    CB1 = 300
    CH1 = 300
    CB2 = 400
    CH2 = 400

    unitVolume1 = (2*CB1+2*CH1)*(CH1)*tf
    unitVolume2 = (2*CB2+2*CH2)*(CH2)*tf

    unitVolume3 = (2*CB1+2*CH1)*Height*tf
    unitVolume4 = (2*CB2+2*CH2)*Height*tf

    FRPlimit = [3,3,2,1,1,1]

    y[0][0] = 0
    for i in range(1,nVar+1):
        if i < 4:
            y[0][0] = y[0][0] + unitVolume1*x[i-1]*2
        else:
            y[0][0] = y[0][0] + unitVolume2*x[i-1]*2
        
        
        if x[i-1] >= FRPlimit[i-1]:
            if i < 4:
                y[0][0] = y[0][0] + ( unitVolume3 - unitVolume1 ) *2
            else:
                y[0][0] = y[0][0] + ( unitVolume4 - unitVolume2 ) *2
            
    Individual = np.zeros((2*nVar,1))
    for i in range(1,nVar+1):
        if i <= nStory:
            Individual[i-1][0] = x[i-1]
            Individual[i+nStory-1][0] = x[i-1]
        else:
            Individual[i+1*nStory-1][0] = x[i-1]
            Individual[i+2*nStory-1][0] = x[i-1]
        

    for i in range(1,nVar*2+1):
        if Individual[i-1][0] == 1:
            Individual[i-1][0] = 0

    nFRP = Individual

    # Fiber definition
    # Steel fiber
    fy = 300.0
    Es = 200000

    hardeningratio = 0.03

    RebarSpacings = 25.0
    ShearSpacings = 300

    # Concrete fiber
    BSec1 = 300
    HSec1 = 300

    BSec2 = 400
    HSec2 = 400

    BSec3 = 250
    HSec3 = 500

    A1 = HSec1*BSec1
    A2 = HSec2*BSec2
    A3 = HSec3*BSec3

    I1 = CRC*HSec1*HSec1*HSec1*BSec1/12.0
    I2 = CRC*HSec2*HSec2*HSec2*BSec2/12.0
    I3 = CRB*HSec3*HSec3*HSec3*BSec3/12.0

    D10 = 9.53
    A10 = 71.3

    D19 = 19.1
    A19 = 286.5

    D22 = 22.2 
    A22 = 387.1

    fc = 21.0
    Ec = 4700*math.sqrt(fc)    
    eps1 = 2*fc/Ec
    eps2 = 0.0035

    # Fiber section
    numP = 6
    numBarsSec11 = 4
    numBarsSec12 = 2
    numBarsSec13 = 2
    numBarsSec14 = 4

    numBarsSec21 = 5
    numBarsSec22 = 2
    numBarsSec23 = 2
    numBarsSec24 = 2
    numBarsSec25 = 5

    numBarsSec31 = 3
    numBarsSec32 = 2

    numBarsSec41 = 3
    numBarsSec42 = 2
    numBarsSec43 = 2

    cover = 40
    coverY1 = HSec1/2.0
    coverZ1 = BSec1/2.0
    coreY1 = coverY1 - cover - D10 - D19/2
    coreZ1 = coverZ1 - cover - D10 - D19/2

    coverY2 = HSec2/2.0
    coverZ2 = BSec2/2.0
    coreY2 = coverY2 - cover - D10 - D22/2
    coreZ2 = coverZ2 - cover - D10 - D22/2

    coverY3 = HSec3/2.0
    coverZ3 = BSec3/2.0
    coreY3 = coverY3 - cover - D10 - D19/2
    coreZ3 = coverZ3 - cover - D10 - D19/2

    nfY = 20
    nfZ = 20

    # Create nodes
    for i in range(1,nStory+2):
        for j in range(1,nSpan+2):
            node(i*10 + j, Span*(j-1), Height*(i-1))

    for j in range(1,nSpan+2):
        i = 1
        fix(i*10 + j, 1, 1, 1)        
        
    for i in range(2,nStory+2):
        for j in range(2,nSpan+2):
            equalDOF(i*10 + 1, i*10+j, 1)

    # mass
    StoryTotalWeight = unitWeight*(nSpan*Span)*EffectiveDepth
    SpanWeight = StoryTotalWeight/nSpan

    for i in range(2,nStory+2):
        for j in range(1,nSpan+2):
            
            if j==1 or j==nSpan+1:
                mass(i*10 + j, SpanWeight[i-2][0]/g/2,0,0)  
            else:
                mass(i*10 + j, SpanWeight[i-2][0]/g,0,0) 

    # Shear definition (assumed very large for reinforced ones)
    G = 0.4*Ec
    Kshearspring01 = (5/6)*BSec1*HSec1*G/Height
    Kshearspring02 = (5/6)*BSec2*HSec2*G/Height
    BigNumber = 10000000000
    nColumn = (nStory)*(nSpan+1)
    Vn = BigNumber

    for i in range(1,7):
        
        if i <= 3:
            j = i
            ShearPoint[i-1][0] = Vn
            ShearPoint[i-1][1] = Vn/Kshearspring01/HSec1
            ShearPoint[i-1][2] = Vn*0.3
            ShearPoint[i-1][3] = Vn/Kshearspring01/HSec1*10
        else:
            j=i+6
            ShearPoint[i-1][0]= Vn
            ShearPoint[i-1][1] = Vn/Kshearspring02/HSec2
            ShearPoint[i-1][2] = Vn*0.3
            ShearPoint[i-1][3] = Vn/Kshearspring02/HSec2*10
            
        uniaxialMaterial('Hysteretic', 500+i, ShearPoint[i-1][0], ShearPoint[i-1][1], ShearPoint[i-1][2], ShearPoint[i-1][3], (-1)*ShearPoint[i-1][0], (-1)*ShearPoint[i-1][1], (-1)*ShearPoint[i-1][2], (-1)*ShearPoint[i-1][3], 0, 0, 0, 0)

    # FRP Effects on existing materials 
    uniaxialMaterial('Steel01',1,fy,Es,hardeningratio)
    uniaxialMaterial('Concrete01',2,-fc,-eps1,-fc,-eps2)

    for i in range(1,nColumnMaterial+1):
        
        niFRP = nFRP[i-1][0]
    #     if i<=(2*nStory):
    #         Ast = ( numBarsSec11 + numBarsSec12 + numBarsSec13 + numBarsSec14 ) * A19
    #     else:        
                
        if nFRP[i-1][0] < 0.1:
            uniaxialMaterial('Concrete01',100+i, fc*(-1), eps1*(-1), fc*(-1), eps2*(-1))
            Eccu[i-1][0] = eps2
        else:
            if (i<=nStory*2):
                # Concrete info
                Ast = ( numBarsSec11 + numBarsSec12 + numBarsSec13 + numBarsSec14 ) * A19
                b = BSec1
                h = HSec1
                rc = 25
                tf = 0.33
                efu_ = 0.0167
                Ef = 227527
                n = niFRP
                Ce = 0.95
                Ke = 0.55
                Pf = 0.95
                ec = eps1
                Ag = h*b
                rhog = Ast/Ag
                D = math.sqrt(b**2 + h**2)
                AreaRatio = (1- ( (b/h)*((h-2*rc)**2) + (h/b)*((b-2*rc)**2) )/(3*Ag) - rhog) / (1-rhog)
                Ka = AreaRatio*((b/h)**2)
                Kb = AreaRatio*((h/b)**0.5)
                efu = Ce*efu_
                efe = Ke*efu
                
                if efe > 0.004:
                    efe = 0.004

                fl = (2*Ef*n*tf*efe) / D
                Min_n = 0.08*fc/(2*Ef*tf*efe)*D
                Ec = 4700*math.sqrt(fc)
                eccu = ec * (1.5 + 12*Kb*(fl/fc)*((efe/ec)**0.45))
                            
                if eccu >= 0.01:
                    eccu = 0.01
                            
                fcc = fc + Pf*3.3*Ka*fl
                E2 = (fcc-fc)/eccu
                et = (2*fc)/(Ec-E2)
                tmpfc = fc + E2*et
                tmpet = et
                tmpfcc = fcc
                tmpeccu = eccu
                uniaxialMaterial('Concrete01',100+i, (-1)*tmpfc, (-1)*tmpet, (-1)*tmpfcc, (-1)*tmpeccu )
            else:
                # Concrete info
                Ast = ( numBarsSec21 + numBarsSec22 + numBarsSec23 + numBarsSec24 + numBarsSec25 ) * A22
                b = BSec2
                h = HSec2
                rc = 25            
                tf = 0.33
                efu_ = 0.0167
                Ef = 227527
                n = niFRP
                Ce = 0.95
                Ke = 0.55
                Pf = 0.95
                ec = eps1
                Ag = h*b
                rhog = Ast/Ag
                D = math.sqrt(b**2 + h**2)
                AreaRatio = (1- ( (b/h)*((h-2*rc)**2) + (h/b)*((b-2*rc)**2) )/(3*Ag) - rhog) / (1-rhog)
                Ka = AreaRatio*((b/h)**2)
                Kb = AreaRatio*((h/b)**0.5)
                efu = Ce*efu_
                efe = Ke*efu
                if efe > 0.004:
                    efe = 0.004

                fl = (2*Ef*n*tf*efe) / D
                Min_n = 0.08*fc/(2*Ef*tf*efe)*D
                Ec = 4700*math.sqrt(fc)
                eccu = ec * (1.5 + 12*Kb*(fl/fc)*((efe/ec)**0.45))
                if eccu > 0.01:
                    eccu = 0.01
                else:
                    eccu = eccu
                            
                fcc = fc + Pf*3.3*Ka*fl
                E2 = (fcc-fc)/eccu
                et = (2*fc)/(Ec-E2)
                tmpfc = fc + E2*et
                tmpet = et
                tmpfcc = fcc
                tmpeccu = eccu
                uniaxialMaterial('Concrete01',100+i, (-1)*tmpfc, (-1)*tmpet, (-1)*tmpfcc, (-1)*tmpeccu )
            
            Eccu[i-1][0] = tmpeccu

    # Insert fiber
    nColumnSection = (nStory)*(nSpan+1)
    tmp1 = (nStory*2)

    for i in range(1,nColumnSection+1):
        if i <= tmp1:
            Ccover = cover + D10 + D19/2
            section('Fiber',i)
            patch('quad',100+i,nfZ,nfY,-HSec1/2, BSec1/2, -HSec1/2, -BSec1/2,HSec1/2, -BSec1/2,HSec1/2, BSec1/2)
            layer('straight', 1, numBarsSec11, A19, HSec1/2-Ccover,       BSec1/2-Ccover, HSec1/2-Ccover,       -(BSec1/2-Ccover))
            layer('straight', 1, numBarsSec12, A19, (HSec1/2-Ccover)*1/3, BSec1/2-Ccover, (HSec1/2-Ccover)*1/3, -(BSec1/2-Ccover))
            layer('straight', 1, numBarsSec13, A19, -(HSec1/2-Ccover)*1/3,BSec1/2-Ccover, -(HSec1/2-Ccover)*1/3,-(BSec1/2-Ccover))
            layer('straight', 1, numBarsSec14, A19, -(HSec1/2-Ccover),    BSec1/2-Ccover, -(HSec1/2-Ccover),    -(BSec1/2-Ccover))        
        #else:
        #    Ccover = cover + D10 + D22/2
        #    section('Fiber',i)
        #    patch('quad',100+i,nfZ,nfY,-HSec2/2, BSec2/2, -HSec2/2, -BSec2/2,HSec2/2, -BSec2/2,HSec2/2, BSec2/2)
        #    layer('straight', 1, numBarsSec21, A22, HSec1/2-Ccover,       BSec1/2-Ccover, HSec1/2-Ccover,       -(BSec1/2-Ccover))
        #    layer('straight', 1, numBarsSec22, A22, (HSec1/2-Ccover)*1/2, BSec1/2-Ccover, (HSec1/2-Ccover)*1/2, -(BSec1/2-Ccover))
        #    layer('straight', 1, numBarsSec23, A22, 0,       BSec1/2-Ccover, 0,       -(BSec1/2-Ccover))
        #    layer('straight', 1, numBarsSec24, A22, -(HSec1/2-Ccover)*1/2,BSec1/2-Ccover, -(HSec1/2-Ccover)*1/2,-(BSec1/2-Ccover))
        #    layer('straight', 1, numBarsSec25, A22, -(HSec1/2-Ccover),    BSec1/2-Ccover, -(HSec1/2-Ccover),    -(BSec1/2-Ccover))
        else:
            Ccover = cover + D10 + D22/2
            section('Fiber',i)
            patch('quad',100+i,nfZ,nfY,-HSec2/2, BSec2/2, -HSec2/2, -BSec2/2,HSec2/2, -BSec2/2,HSec2/2, BSec2/2)
            layer('straight', 1, numBarsSec21, A22, HSec2/2-Ccover,       BSec2/2-Ccover, HSec2/2-Ccover,       -(BSec2/2-Ccover))
            layer('straight', 1, numBarsSec22, A22, (HSec2/2-Ccover)*1/2, BSec2/2-Ccover, (HSec2/2-Ccover)*1/2, -(BSec2/2-Ccover))
            layer('straight', 1, numBarsSec23, A22, 0,       BSec2/2-Ccover, 0,       -(BSec2/2-Ccover))
            layer('straight', 1, numBarsSec24, A22, -(HSec2/2-Ccover)*1/2,BSec2/2-Ccover, -(HSec2/2-Ccover)*1/2,-(BSec2/2-Ccover))
            layer('straight', 1, numBarsSec25, A22, -(HSec2/2-Ccover),    BSec2/2-Ccover, -(HSec2/2-Ccover),    -(BSec2/2-Ccover))    


    ExtShearSpringData = [1,2,3,1,2,3]
    IntShearSpringData = [6,4,5]
    for i in range(1,nColumnSection+1):
        
        if i <= 6:
            #section('Aggregator',100+i, 500+ExtShearSpringData[i-1],'Vy','-section', i)
            section('Aggregator',100+i, 501,'Vy','-section', i)
        else:
            j= (i % 3)+1
            #section('Aggregator',100+i, 500+IntShearSpringData[j-1],'Vy','-section', i)
            section('Aggregator',100+i, 501,'Vy','-section', i)

    section('Fiber',53)
    patch('quad', 2,nfZ,nfY,-coverY3, coverZ3, -coverY3, -coverZ3, coverY3, -coverZ3, coverY3, coverZ3)
    layer('straight', 1,numBarsSec31,A19, coreY3, coreZ3,  coreY3, -coreZ3)
    layer('straight', 1,numBarsSec32,A19,-coreY3, coreZ3, -coreY3, -coreZ3)

    section('Fiber', 54)
    patch('quad', 2, nfZ, nfY, -coverY3, coverZ3, -coverY3, -coverZ3, coverY3, -coverZ3, coverY3, coverZ3)
    layer('straight', 1, numBarsSec41, A19,  coreY3,coreZ3,  coreY3, -coreZ3)
    layer('straight', 1, numBarsSec42, A19,  coreY3-25-D19,  coreZ3, coreY3-25-D19, -coreZ3)
    layer('straight', 1, numBarsSec43, A19, -coreY3, coreZ3, -coreY3, -coreZ3)


    # Defining beam with hinges
    nColumn = (nStory)*(nSpan+1)
    geomTransf('Linear',2)
    geomTransf('PDelta',1)
    ##section('Fiber',501)
    section('Elastic',990501,Ec,A1,I1)
    ##section('Fiber',502)
    section('Elastic',990502,Ec,A2,I2)
    section('Elastic',990503,Ec,A3,I3)

    for i in range(1,nColumn+1):
        iStory = (i % nStory)
        if iStory == 0:
            iStory = 3
            
        iSpan = math.floor( (i-1)/nStory )
        
        if i <= nStory:
            #element('beamWithHinges',i,10*iStory+1,10*(iStory+1)+1,100+i,HSec1*0.5,100+i,HSec1*0.5,Ec,A1,I1, 1)
            beamIntegration('HingeRadau',10000+i, 100+i,HSec1*0.5,100+i,HSec1*0.5,990501)
            element('forceBeamColumn',i,10*iStory+1,10*(iStory+1)+1, 1,10000+i)
        elif i <= (2*nStory):
            #element('beamWithHinges', i, 10*iStory+(nSpan+1), 10*(iStory+1)+(nSpan+1), 100+i, HSec1*0.5, 100+i, HSec1*0.5,Ec,A1,I1,1)
            beamIntegration('HingeRadau',10000+i, 100+i,HSec1*0.5,100+i,HSec1*0.5,990501)
            element('forceBeamColumn',i,10*iStory+(nSpan+1), 10*(iStory+1)+(nSpan+1), 1,10000+i)
        else:
            #element('beamWithHinges', i, 10*iStory+iSpan, 10*(iStory+1)+iSpan, 100+i,HSec2*0.5, 100+i, HSec2*0.5,Ec,A2,I2,1)
            beamIntegration('HingeRadau',10000+i, 100+i,HSec1*0.5,100+i,HSec1*0.5,990502)
            element('forceBeamColumn',i,10*iStory+iSpan, 10*(iStory+1)+iSpan, 1,10000+i)     


    # Manual definition
    beamIntegration('HingeRadau',10000+13, 53, HSec3*0.5,54,HSec3*0.5,990503)
    element('forceBeamColumn',13,21,22,2,10000+13)
    beamIntegration('HingeRadau',10000+14, 54, HSec3*0.5,54,HSec3*0.5,990503)
    element('forceBeamColumn',14,22,23,2,10000+14)
    beamIntegration('HingeRadau',10000+15, 54, HSec3*0.5,53,HSec3*0.5,990503)
    element('forceBeamColumn',15,23,24,2,10000+15)

    beamIntegration('HingeRadau',10000+16, 53, HSec3*0.5,54,HSec3*0.5,990503)
    element('forceBeamColumn',16,31,32,2,10000+16)
    beamIntegration('HingeRadau',10000+17, 54, HSec3*0.5,54,HSec3*0.5,990503)
    element('forceBeamColumn',17,32,33,2,10000+17)
    beamIntegration('HingeRadau',10000+18, 54, HSec3*0.5,53,HSec3*0.5,990503)
    element('forceBeamColumn',18,33,34,2,10000+18)

    beamIntegration('HingeRadau',10000+19, 53, HSec3*0.5,54,HSec3*0.5,990503)
    element('forceBeamColumn',19,41,42,2,10000+19)
    beamIntegration('HingeRadau',10000+20, 54, HSec3*0.5,54,HSec3*0.5,990503)
    element('forceBeamColumn',20,42,43,2,10000+20)
    beamIntegration('HingeRadau',10000+21, 54, HSec3*0.5,53,HSec3*0.5,990503)
    element('forceBeamColumn',21,43,44,2,10000+21)


    # Analysis
    timeSeries("Linear",2)
    pattern('Plain',1,2)
    eleLoad('-ele', 13, 14, 15, '-type', '-beamUniform', beamload)
    eleLoad('-ele', 16, 17, 18, '-type', '-beamUniform', beamload)
    eleLoad('-ele', 19, 20, 21, '-type', '-beamUniform', beamload)

    #DisplayModel2D NodeNumbers
    #DisplayModel2D DeformedShape $ViewScale ;	# display deformed shape, the scaling factor needs to be adjusted for each model

    numEigen = 1  # Only the first mode
    eigenValues = eigen(numEigen)
    omega1 = eigenValues[0]**0.5
    T1 = 2 * 3.14159 / omega1

    system('BandGeneral')
    constraints('Plain')
    numberer('RCM')
    test('EnergyIncr', 1.0e-8, 50)
    algorithm('Newton')
    integrator('LoadControl', 0.1, 5, 0.1, 0.1)
    analysis('Static')
    analyze(10)

    #print("Gravity load analysis completed")

    loadConst('-time', 0.0)

    pattern('Plain',3,2)
    load(21, 0.167, 0.0, 0.0)
    load(31, 0.333, 0.0, 0.0)
    load(41, 0.5, 0.0, 0.0)

    #recorder('Drift', '-file', "results/EQDrift.out", '-time', '-iNode', 11, 21, 31, '-jNode', 21, 31, 41, '-dof', 1, '-perDirn', 2)
    recorder('Node', '-file', "results/node.out", '-time', '-node', 41, '-dof', 1, 'disp')
    recorder('Node', '-file', "results/nodeDrift.out", '-time', '-node', 11, 21, 31, 41, '-dof', 1, 'disp')
    recorder('Node', '-file', "results/Vbase.out", '-node', 11, 12, 13, 14, '-dof', 1, 'reaction')

    recorder('Element', '-file', "results/ColumnEleSection6force.out", '-time', '-ele', 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 'section',  6, 'force')
    recorder('Element', '-file', "results/ColumnEleSection6Deformation.out", '-time', '-ele', 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 'section', 6, 'deformation')
    recorder('Element', '-file', "results/ColumnEleSection1force.out", '-time', '-ele', 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 'section',  1, 'force')
    recorder('Element', '-file', "results/ColumnEleSection1Deformation.out", '-time', '-ele', 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 'section',  1, 'deformation')

    recorder('Element', '-file', "results/BeamEleSection6force.out", '-time', '-ele', 13, 14, 15, 16, 17, 18, 19, 20, 21, 'section',  6, 'force')
    recorder('Element', '-file', "results/BeamEleSection6Deformation.out", '-time', '-ele', 13, 14, 15, 16, 17, 18, 19, 20, 21, 'section',  6, 'deformation')
    recorder('Element', '-file', "results/BeamEleSection1force.out", '-time', '-ele', 13, 14, 15, 16, 17, 18, 19, 20, 21, 'section',  1, 'force')
    recorder('Element', '-file', "results/BeamEleSection1Deformation.out", '-time', '-ele', 13, 14, 15, 16, 17, 18, 19, 20, 21, 'section',  1, 'deformation')

    recorder('Element','-file', "results/EleLocalforce.out", '-time', '-ele', 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 'localForce')

    recorder('Element', '-file', "results/Strain1.out", '-time', '-ele', 1, 2, 3, 4, 5, 6, 'section', 1, 'fiber',coverY1, 0, 'stressStrain')
    recorder('Element', '-file', "results/Strain2.out", '-time', '-ele', 1, 2, 3, 4, 5, 6, 'section', 1, 'fiber', -coverY1, 0, 'stressStrain')
    recorder('Element', '-file', "results/Strain3.out", '-time', '-ele', 1, 2, 3, 4, 5, 6, 'section', numP, 'fiber', coverY1, 0, 'stressStrain')
    recorder('Element', '-file', "results/Strain4.out", '-time', '-ele', 1, 2, 3, 4, 5, 6, 'section', numP, 'fiber', -coverY1, 0, 'stressStrain')

    recorder('Element', '-file', "results/Strain5.out", '-time', '-ele', 7, 8, 9, 10, 11, 12, 'section', 1, 'fiber', coverY2, 0, 'stressStrain')
    recorder('Element', '-file', "results/Strain6.out", '-time', '-ele', 7, 8, 9, 10, 11, 12, 'section', 1, 'fiber', -coverY2, 0, 'stressStrain')
    recorder('Element', '-file', "results/Strain7.out", '-time', '-ele', 7, 8, 9, 10, 11, 12, 'section', numP, 'fiber', coverY2, 0, 'stressStrain')
    recorder('Element', '-file', "results/Strain8.out", '-time', '-ele', 7, 8, 9, 10, 11, 12, 'section', numP, 'fiber', -coverY2, 0, 'stressStrain')


    ## constraintsType = Plain;      # default; 
    # constraints('Plain')
    # numberer('RCM')
    ## system Umfpack
    # system('BandGeneral') 

    Tol = 1.e-6        # Convergence Test: tolerance 
    maxNumIter = 20      # Convergence Test: maximum number of iterations that will be performed before "failure to converge" is returned 
    printFlag  = 0      # Convergence Test: flag used to print information on convergence (optional)        # 1: print information on each step; 
    #set TestType NormDispIncr;   # Convergence-test type 
    test('NormDispIncr', Tol, maxNumIter, printFlag) 

    #set algorithmType Newton 
    #set algorithm KrylovNewton
    algorithm('Newton')

    Dincr_min  = Dincr*1.0
    Dincr_max  = Dincr*1.0
    Jd = 4

    integrator('DisplacementControl', IDctrlNode,IDctrlDOF,Dincr,Jd,Dincr_min,Dincr_max)
    analysis('Static')

    #  ---------------------------------    perform Static Pushover Analysis 
    Nsteps = int(Dmax/Dincr)        # number of pushover analysis steps 
    ok = analyze(Nsteps)                # this will return zero if no convergence problems were encountered 

    if ok != 0:      
    # if analysis fails, we try some other stuff, performance is slower inside this loop 
        ok  = 0     
        controlDisp = 0.0
        D0 = 0.0      # analysis starts from zero
        Dstep = (controlDisp-D0)/(Dmax-D0) 

        while Dstep < 1.0 and ok == 0:    
        
            controlDisp = nodeDisp(IDctrlNode,IDctrlDOF) 
            Dstep = (controlDisp-D0)/(Dmax-D0) 
            ok = analyze(1) 
            if ok != 0: 
                print("Trying Newton with Initial Tangent ..") 
                test('NormDispIncr',Tol, 2000,  0) 
                algorithm('Newton','-initial') 
                ok = analyze(1) 
                test('NormDispIncr',Tol,maxNumIter,  0) 
                algorithm('Newton') 
        
            if ok != 0: 
                print("Trying Broyden ..") 
                algorithm('Broyden', 8) 
                ok = analyze(1) 
                algorithm('Newton') 
        
            if ok != 0: 
                print("Trying NewtonWithLineSearch ..") 
                algorithm('NewtonLineSearch', .8) 
                ok = analyze(1) 
                algorithm('Newton')

    #print("Pushover analysis completed!!!")

    TotalWeight = sum(StoryTotalWeight)
    [DispIO, DispLS, DispCP, Area02] = TargetDisp(TotalWeight,T1)

    if DispIO > nStory*Height*0.01:
        DispIO = nStory*Height*0.01*0.98
        
    if DispLS > nStory*Height*0.02:
        DispLS = nStory*Height*0.02*0.98
        
    if DispCP > nStory*Height*0.04:
        DispCP = nStory*Height*0.04*0.98

    nodeIO, nodeLS, nodeCP, DriftIO, DriftLS, DriftCP, IndexStep = ReadOutput(DispIO, DispLS, DispCP, nStory) # check ReadOutput.m
    tmpCons = DriftLS[0]/AcceptanceCriteria[1]

    # objective function 1
    y[0][0] = y[0][0] / 30492000
    
    # objective function 2
    y[0][1] = 1/Area02*10**9
    
    if nObj == 3:
        # objective function 3
        individual3 = np.zeros(2*nVar)
        N_re = 0
        Story = np.array([1,2,3,1,2,3])
        for i in range(nVar):
            if i < nStory:
                individual3[i] = x[i]
                individual3[i+nStory] = x[i]
                N_re = N_re + Story[i]*individual3[i]*2
            else:
                individual3[i+1*nStory] = x[i]
                individual3[i+2*nStory] = x[i]
                N_re = N_re + Story[i] * individual3[i+1*nStory]*2
        y[0][2] = N_re/(6*8*6)*1.0

        # constraints
        # ConsValue = EvalConstraint(nCon, nVar, IndexStep, Eccu) # check EvalConstraint.m
        # cons1 = np.append([tmpCons],ConsValue)
        ConsValue = EvalConstraint(nCon+1, nVar, IndexStep, Eccu) # check EvalConstraint.m
        cons1 = ConsValue
        
    elif nObj == 4:
        # objective function 3
        individual3 = np.zeros(2*nVar)
        N_re = 0
        Story = np.array([1,2,3,1,2,3])
        for i in range(nVar):
            if i < nStory:
                individual3[i] = x[i]
                individual3[i+nStory] = x[i]
                N_re = N_re + Story[i]*individual3[i]*2
            else:
                individual3[i+1*nStory] = x[i]
                individual3[i+2*nStory] = x[i]
                N_re = N_re + Story[i] * individual3[i+1*nStory]*2
        y[0][2] = N_re/(6*8*6)*1.0
        
        # objective function 4
        y[0][3] = tmpCons
        
        # constraints
        ConsValue = EvalConstraint(nCon+1, nVar, IndexStep, Eccu) # check EvalConstraint.m
        cons1 = ConsValue
    
    for i in range(len(cons1)):
        if cons1[i] > 1.0:
            cons1[i] = cons1[i] - 1.0
        else:
            cons1[i] = 0.0
    
    cons =  np.round(cons1*100)/100

    return y, cons

'''
Utility functions for computing objectives and constraints
'''

def EvalConstraint(nCon, nVar, IndexStep, Eccu):
    nVar = nVar*2
    nCon1 = nCon-1
    ConsValue = np.zeros(nCon1)
    tmpConsValue = np.zeros(nCon1)
    nStep = int(IndexStep[1]) # LS

    ## strain (concrete fiber)
    Strain1 = pd.read_csv("results/Strain1.out", sep=" ", header=None).to_numpy()
    Strain2 = pd.read_csv("results/Strain2.out", sep=" ", header=None).to_numpy()
    Strain3 = pd.read_csv("results/Strain3.out", sep=" ", header=None).to_numpy()
    Strain4 = pd.read_csv("results/Strain4.out", sep=" ", header=None).to_numpy()

    Strain5 = pd.read_csv("results/Strain5.out", sep=" ", header=None).to_numpy()
    Strain6 = pd.read_csv("results/Strain6.out", sep=" ", header=None).to_numpy()
    Strain7 = pd.read_csv("results/Strain7.out", sep=" ", header=None).to_numpy()
    Strain8 = pd.read_csv("results/Strain8.out", sep=" ", header=None).to_numpy()

    maxE = np.zeros(nVar)
    StrainValue = np.zeros(nVar)

    for i in range(nVar):
        if i < nVar/2.0 :
            j=i
            tmp1 = Strain1[nStep, 2*(j+1)]
            tmp2 = Strain2[nStep, 2*(j+1)]
            tmp3 = Strain3[nStep, 2*(j+1)]
            tmp4 = Strain4[nStep, 2*(j+1)]
        else:
            j=i-2*3
            tmp1 = Strain5[nStep, 2*(j+1)]
            tmp2 = Strain6[nStep, 2*(j+1)]
            tmp3 = Strain7[nStep, 2*(j+1)]
            tmp4 = Strain8[nStep, 2*(j+1)]

        maxE[i] = np.min([tmp1,tmp2,tmp3,tmp4])
        StrainValue[i] = np.abs( maxE[i]/Eccu[i].item() )
        tmpConsValue[i] = StrainValue[i]
        ConsValue[i] = tmpConsValue[i]

    return ConsValue

def ReadOutput(DispIO, DispLS, DispCP, nStory):

    node = pd.read_csv("results/node.out", sep=" ", header=None).to_numpy()
    nodeDrift = pd.read_csv("results/nodeDrift.out", sep=" ", header=None).to_numpy()

    m, n = nodeDrift.shape
    
    EQDrift = np.zeros((m,n-1))
    EQDrift[:,0] = nodeDrift[:,0]
    EQDrift[:,1:] = np.diff(nodeDrift[:,1:])/3000
    
    m, n = node.shape

    IndexStep = np.zeros(3)
    
    for i in range(1,m):
        if ( (node[i-1,1]<=DispIO) & (node[i,1]>=DispIO) ):
            IndexStep[0] = i
            nodeIO = node[i,:]
            DriftIO = EQDrift[i,:]
            DriftIO[0] = max(EQDrift[i,1:])
    
        if ( (node[i-1,1]<=DispLS) & (node[i,1]>=DispLS) ):  
            IndexStep[1] = i
            nodeLS = node[i,:]
            DriftLS = EQDrift[i,:]
            DriftLS[0] = max(EQDrift[i,1:])    
    
        if ( (node[i-1,1]<=DispCP) & (node[i,1]>=DispCP) ):  
            IndexStep[2] = i
            nodeCP = node[i,:]
            DriftCP = EQDrift[i,:]
            DriftCP[0] = max(EQDrift[i,1:])  
            tmp = DriftCP*1000
            DriftCP = np.round(tmp)/1000
            break

    return nodeIO, nodeLS, nodeCP, DriftIO, DriftLS, DriftCP, IndexStep


def TargetDisp(TotalWeight,Ti):

    # 1. C0
    C0 = 1.2 # Story 3

    # Ki : Elastic transverse stiffness of the building
    # Ke : Effective lateral stiffness of the building
    Ke, Ki, YieldStrength, postelasticstiffness, Area02 = Bi_Linear2()
    #load Period.out -ascii

    #Ti=[Period]
    
    # effective period of the building
    Te = Ti*np.sqrt(Ki/Ke)

    # C1
    Ss = np.array([0.4121,0.9445,1.6220])  # 50%/50, 10%/50, 2%/50
    S1 = np.array([0.1513,0.3468,0.6244])

    Fa = np.array([1.4703,1.1222,1.0])
    Fv = np.array([2.1947,1.7064,1.5])

    Sxs = Fa * Ss
    Sx1 = Fv * S1
    Bs = 1.0 # damping 5%
    B1 = 1.0 # damping 5%

    Ts = (Sx1*Bs)/(Sxs*B1)
    To = 0.2*Ts

    # Sa 
    Sa = np.zeros(3)
    for i in range(3):
        if Te <= To[i]:
            Sa[i] = Sxs[i]*((5/Bs-2)*(Te/Ts[i])+0.4)
        elif Te <= Ts[i]:
            Sa[i] = Sxs[i]/Bs
        elif Te > Ts[i]:
            Sa[i] = Sx1[i]/(B1*Te)

    w = TotalWeight
    Cm = 0.9
    R = (Sa/(YieldStrength/w))*(Cm)

    C1 = np.zeros(3)
    for i in range(3):
        if Te >= Ts[i]:
            C1[i] = 1.0
        else:
            C1[i]=( 1.0 + (R[i]-1)*Ts[i]/Te ) / R[i]
        
        if Te < 0.1:
            C11 = 1.5

        if (Te >= 0.1) & (Te < Ts[i]):
            C11 = ( (1.0-1.5)/(Ts[i]-0.1) ) * (Te-0.1) + 1.5

        if Te >= Ts[i]:
            C11 = 1.0

        if C1[i] > C11:
            C1[i] = C11
        
        if C1[i] < 1.0:
            C1[i]=1.0

    # C2
    C2 = np.zeros(3)
    for i in range(3):
        if i == 0:
            if Te < 0.1:
                C2[i] = 1.0
            elif Te < Ts[i]:
                C2[i] = (1.0 - 1.0)/(Ts[i]-0.1) * (Te-Ts[i]) + 1.0
            else:
                C2[i] = 1.0
        
        elif i == 1:
            if Te < 0.1:
                C2[i] = 1.3
            elif Te < Ts[i]:
                C2[i] = (1.1 - 1.3)/(Ts[i]-0.1) * (Te-Ts[i]) + 1.1
            else:
                C2[i] = 1.1

        else:
            if Te < 0.1:
                C2[i] = 1.5
            elif Te < Ts[i]:
                C2[i] = (1.2 - 1.5)/(Ts[i]-0.1) * (Te-Ts[i]) + 1.2
            else:
                C2[i] = 1.2
        
        if i == 1:
            C2[i] = 1.0

    # C3
    stiffnessratio=postelasticstiffness/Ke
    C3 = np.zeros(3)
    for i in range(3):
        if postelasticstiffness < 0:
            C3[i]=1.0 + abs(stiffnessratio)*(R[i]-1)^(3/2)/Te
        else:
            C3[i]=1.0

    g=9806.65 #mm/sec^2

    DispIO = C0*C1[0]*C2[0]*C3[0]*Sa[0]*(Te**2)*g/(4*math.pi**2)
    DispLS = C0*C1[1]*C2[1]*C3[1]*Sa[1]*(Te**2)*g/(4*math.pi**2)
    DispCP = C0*C1[2]*C2[2]*C3[2]*Sa[2]*(Te**2)*g/(4*math.pi**2)

    return DispIO, DispLS, DispCP, Area02

def Bi_Linear2():

    node = pd.read_csv("results/node.out", sep=" ", header=None).to_numpy()
    Pushover = np.vstack(([0,0],node))
    m, n = Pushover.shape
    Vmax = max(Pushover[:,0])
    BiLinear = np.zeros((3,2))
    K = 0.0

    # Area from original pushover curve
    Area01 = 0.0
    for i in range(m-1):
        Area01 += (Pushover[i+1,1]-Pushover[i,1])*(Pushover[i+1,0] + Pushover[i,0])/2.0

    Vy = Vmax
    flag01 = True
    while flag01:
        Error = 0

        Vy06 = 0.6*Vy
        i = 0
        flag02 = True
        while flag02:
            i += 1
            if Pushover[i,0] >= Vy06:
                flag02 = False
        
                K = Pushover[i,0]/Pushover[i,1]
                BiLinear[1,0] = Vy/K
                BiLinear[1,1] = Vy
                BiLinear[2,0] = Pushover[m-1,1]
                BiLinear[2,1] = Pushover[m-1,0]
            
                Area02 = 0.5*BiLinear[1,0]*BiLinear[1,1] + (0.5*(BiLinear[2,0]-BiLinear[1,0])*(BiLinear[2,1]+BiLinear[1,1]))
        
        if abs(Area01-Area02) <= (0.01*Area01):
            flag01 = False

        if Vy <= 0.6*Vmax:
            Error = 1
            flag01 = False

        Vy -= Vmax*0.01

    Ke = K
    Ki = Pushover[2,0]/Pushover[2,1]
    YieldStrength = Vy
    postelasticstiffness = (BiLinear[2,1]-BiLinear[1,1])/(BiLinear[2,0]-BiLinear[1,0])

    return Ke, Ki, YieldStrength, postelasticstiffness, Area02
