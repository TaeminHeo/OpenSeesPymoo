function [DispIO, DispLS, DispCP,Area02]=TargetDisp(TotalWeight)

% 목표변위 

%1.C0
C0= 1.2; %3층일때 

% Ki : 건물의 탄성횡강성 Ke :건물의 유효 횡강성
[Ke,Ki,YieldStrength,postelasticstiffness,Area02]=Bi_Linear2();
% 탄성동적해석으로 계산된 탄성 1차 모드 주기
% load Period.out -ascii; 

Ti=[3.06605237232262];

% effective period of the building
Te=Ti*(Ki/Ke)^(1/2);

% C1
Ss = [0.4121; 0.9445; 1.6220];  % 50%/50, 10%/50, 2%/50
S1 = [0.1513; 0.3468; 0.6244];

Fa = [1.4703; 1.1222; 1.0];
Fv = [2.1947; 1.7064; 1.5];

Sxs=Fa.*Ss;
Sx1=Fv.*S1;
Bs=1.0; % damping 5%
B1=1.0; % damping 5%

Ts = (Sx1*Bs)./(Sxs*B1);
To = 0.2*Ts;

% Sa 
Sa = zeros(3,1);
for i=1:3
    if Te<=To(i)
        Sa(i) = Sxs(i)*((5/Bs-2)*(Te/Ts(i))+0.4);
    elseif Te<=Ts(i)
        Sa(i) = Sxs(i)/Bs;
    elseif Te>Ts(i)
        Sa(i) = Sx1(i)/(B1*Te);
    end
end

w=TotalWeight;
Cm = 0.9;
R= (Sa./(YieldStrength/w))*(Cm);

for i=1:3
    if Te >= Ts(i)
        C1(i) = 1.0;
    else
        C1(i)=( 1.0 + (R(i)-1)*Ts(i)/Te ) / R(i);
    end
     
    if Te<0.1
        C11=1.5;
    end

    if Te>=0.1 && Te<Ts(i)
        C11 = ( (1.0-1.5)/(Ts(i)-0.1) ) * (Te-0.1) + 1.5;
    end

    if Te>=Ts(i)
        C11=1.0;
    end

    if C1(i)>C11
        C1(i)=C11;
    end
    
    if C1(i)<1.0
        C1(i)=1.0;
    end
end

% C2
for i=1:3
    if i == 1
        if Te < 0.1
            C2(i) = 1.0;
        elseif Te < Ts(i)
            C2(i) = (1.0 - 1.0)/(Ts(i)-0.1) * (Te-Ts(i)) + 1.0;
        else
            C2(i) = 1.0;
        end
    elseif i ==2
        if Te < 0.1
            C2(i) = 1.3;
        elseif Te < Ts(i)
            C2(i) = (1.1 - 1.3)/(Ts(i)-0.1) * (Te-Ts(i)) + 1.1;
        else
            C2(i) = 1.1;
        end
    else
        if Te < 0.1
            C2(i) = 1.5;
        elseif Te < Ts(i)
            C2(i) = (1.2 - 1.5)/(Ts(i)-0.1) * (Te-Ts(i)) + 1.2;
        else
            C2(i) = 1.2;
        end
    end
    
    if i == 2
        C2(i) = 1.0;
    end
end

% C3
stiffnessratio=postelasticstiffness/Ke;
for i=1:3
    if postelasticstiffness<0
        C3(i)=1.0 + abs(stiffnessratio)*(R(i)-1)^(3/2)/Te;
    else
        C3(i)=1.0;
    end
end

g=9806.65; %mm/sec^2

DispIO = C0*C1(1)*C2(1)*C3(1)*Sa(1)*(Te^2)*g/(4*pi^2);
DispLS = C0*C1(2)*C2(2)*C3(2)*Sa(2)*(Te^2)*g/(4*pi^2);
DispCP = C0*C1(3)*C2(3)*C3(3)*Sa(3)*(Te^2)*g/(4*pi^2);

