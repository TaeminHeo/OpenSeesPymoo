% from Pushover curve to Bi-linear curve (FEMA 273 based)

function [Ke, Ki,YieldStrength,postelasticstiffness, Area02]=Bi_Linear2()

load('node.out');
Pushover=[0,0;node];
[m,n]=size(Pushover);
[Vmax, Vmax_Index]=max(Pushover(:,1));
BiLinear=zeros(3,2);
K=[0,0];
Pushover(1,1);

% Area from original pushover curve
Area01 = 0.0;               
for i=1:m-1
    Area01 = Area01 + (Pushover(i+1,2)-Pushover(i,2))*(Pushover(i+1,1) + Pushover(i,1))/2.0;
end

Vy=Vmax;
flag01=1;
while (flag01==1)
    Error=0;
    
    Vy06= 0.6*Vy;
    i=1;
    flag02=1;
    while (flag02==1)
       i=i+1;
       if(Pushover(i,1) >= Vy06)
           flag02=0;
       
       
           K(1,1)=Pushover(i,1)/Pushover(i,2);
           BiLinear(2,1)= Vy/K(1,1);
           BiLinear(2,2)= Vy;
           BiLinear(3,1)= Pushover(m,2);
           BiLinear(3,2)= Pushover(m,1);
           
           Area02= 0.5*BiLinear(2,1)*BiLinear(2,2) + ( 0.5*(BiLinear(3,1)-BiLinear(2,1) )*( BiLinear(3,2)+BiLinear(2,2)) );
       end
    end

    if(abs(Area01-Area02) <=(0.01*Area01) )
       flag01=0; 
    end

    if(Vy <= 0.6*Vmax)
        Error=1;
        flag01=0;
    
    end
   
   Vy = Vy- Vmax*0.01;
    
end

Ke=K(1,1);
Ki = Pushover(3,1)/Pushover(3,2);
YieldStrength=Vy;
postelasticstiffness= (BiLinear(3,2)-BiLinear(2,2))/(BiLinear(3,1)-BiLinear(2,1));