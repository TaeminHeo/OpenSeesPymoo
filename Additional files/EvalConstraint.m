function ConsValue = EvalConstraint(nCon, nVar, X, IndexStep, Eccu)
    nVar = nVar*2;
    nCon1 = nCon-1;
    ConsValue = zeros(1,nCon1);
    tmpConsValue = zeros(1,nCon1);
    
    %% calculation: shear strength
    tf = 0.333;
    Ce = 0.95;
    efu_ = 0.0167;
    efu = Ce*efu_;
    efe = 0.75*efu;
    if efe > 0.004
        efe = 0.004;
    end
    
    Ef = 227527;
    
    Vc = zeros(1,nVar);
    Vs = zeros(1,nVar);
    Vf = zeros(1,nVar);
    Vn = zeros(1,nVar);
    
    Nu = [362340	241560	120780	362340	241560	120780	724680	483120	241560	724680	483120	241560];
    efe = 0.004;
    factor1 = 0.95;
    
    fy = 300.0;   % N/mm^2
    
    D10 = 9.53;
    A10 = 71.3;

    D19 = 19.1;
    A19 = 286.5;

    D22 = 22.2; 
    A22 = 387.1;

    BSec1 = 300;  % external column
    HSec1 = 300;

    BSec2 = 400;  % internal column
    HSec2 = 400;

    cover = 40;

    fc = 21.0;
    Ec = 4700*sqrt(fc); 
    
    dd1 = HSec1 - cover - D10 - D19/2;
    dd2 = HSec2 - cover - D10 - D22/2;
    d1 = HSec1;
    d2 = HSec2;
    bw1 = BSec1;
    bw2 = BSec2;
    Ag1 = BSec1*HSec1;
    Ag2 = BSec2*HSec2;
    
    space = 300;
    
    nStep = IndexStep(1,2); % LS
    
    %% strain (concrete fiber)
    load('Strain1.out');
    load('Strain2.out');
    load('Strain3.out');
    load('Strain4.out');
        
    load('Strain5.out');
    load('Strain6.out');
    load('Strain7.out');
    load('Strain8.out');
    
    maxE = zeros(1,nVar);
    StrainValue = zeros(1,nVar);
    StrainValue02 = zeros(nVar,3);
    for i=1:nVar
        if i <= nVar/2.0
            j=i;
            tmp1 = ( Strain1(nStep, 2*j+1) );
            tmp2 = ( Strain2(nStep, 2*j+1) );
            tmp3 = ( Strain3(nStep, 2*j+1) );
            tmp4 = ( Strain4(nStep, 2*j+1) );
        else
            j=i-2*3;
            tmp1 = ( Strain5(nStep, 2*j+1) );
            tmp2 = ( Strain6(nStep, 2*j+1) );
            tmp3 = ( Strain7(nStep, 2*j+1) );
            tmp4 = ( Strain8(nStep, 2*j+1) );
        end
        maxE(i) = ( min([tmp1 tmp2 tmp3 tmp4]) );
        StrainValue(i) = abs( maxE(i)/Eccu(i) );
        StrainValue02(i,:) = [ maxE(i), Eccu(i), abs(maxE(i)/Eccu(i)) ];
        tmpConsValue(i) = StrainValue(i);
        ConsValue(i) = tmpConsValue(i);
    end

    StrainValue02;
end