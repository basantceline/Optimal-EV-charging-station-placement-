clc
clear all
sample_number=1;
w=2*pi*50;


for jj=1:sample_number

% line=[3 0 0 0 0 0
%       0 5 0 0 0 0
%       0 0 6 0 8 0
%       0 0 0 4 0 0
%       0 0 0 0 0 0
%       0 0 0 0 0 7];
% 
% load_power=[0 7 0 2 0 8 1]; %P+jQ

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% line=[1 0 0 0 0 0 0 0
%       0 5 0 0 0 0 4 0
%       0 0 6 0 8 0 0 0
%       0 0 0 4 0 0 0 0
%       0 0 0 0 0 0 0 0
%       0 0 0 0 0 7 0 0
%       0 0 0 0 0 0 0 0
%       0 0 0 0 0 0 0 4];
% 
% load_power=[0 7+j 0 2 3 0 4 0 2];
% source=    [0 0 0 3 9 0 4 0 5];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


line=[(.195+1j*.08)*.6           0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0
            0            (.195+1j*.08)*.55             0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0
            0                    0             (.299+1j*.083)*.55            0                     0                     0                     0                     0                     0                     0                     0              (.524+1j*.090)*.3            0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0      
            0                    0                     0             (.299+1j*.083)*.5             0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                
            0                    0                     0                     0             (.299+1j*.083)*.5             0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0
            0                    0                     0                     0                     0             (.524+1j*.090)*.6             0                     0                     0                     0                     0                     0                     0                     0                     0              (.299+1j*.083)*.6            0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                    
            0                    0                     0                     0                     0                     0             (.524+1j*.090)*.4             0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0              (.524+1j*.090)*.2            0                     0                     0                     0                     0                     0
            0                    0                     0                     0                     0                     0                     0             (.524+1j*.090)*.6             0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0
            0                    0                     0                     0                     0                     0                     0                     0             (.524+1j*.090)*.4             0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0
            0                    0                     0                     0                     0                     0                     0                     0                     0             (.524+1j*.090)*.25            0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0              (.524+1j*.090)*.3            0                     0                     0
            0                    0                     0                     0                     0                     0                     0                     0                     0                     0             (.524+1j*.090)*.2             0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0
            0                    0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0
            0                    0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0              (.524+1j*.090)*.4            0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0
            0                    0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0              (.524+1j*.090)*.2            0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0
            0                    0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0              (.524+1j*.090)*.1            0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0
            0                    0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0
            0                    0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0              (.299+1j*.083)*.55           0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0
            0                    0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0              (.378+1j*.086)*.55           0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0
            0                    0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0              (.378+1j*.086)*.5            0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0
            0                    0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0              (.378+1j*.086)*.5            0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0
            0                    0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0              (.524+1j*.090)*.5            0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0
            0                    0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0              (.524+1j*.090)*.5            0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0
            0                    0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0              (.524+1j*.090)*.6            0                     0                     0                     0                     0                     0                     0                     0                     0                     0
            0                    0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0              (.524+1j*.090)*.4            0                     0                     0                     0                     0                     0                     0                     0                     0
            0                    0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0              (.524+1j*.090)*.25           0                     0                     0                     0                     0                     0                     0                     0
            0                    0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0              (.524+1j*.090)*.2            0                     0                     0                     0                     0                     0                     0
            0                    0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0
            0                    0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0              (.524+1j*.090)*.2            0                     0                     0                     0                     0
            0                    0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0              (.524+1j*.090)*.3            0                     0                     0                     0
            0                    0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0
            0                    0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0             (.524+1j*.090)*.4             0                     0
            0                    0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0             (.524+1j*.090)*.3             0
            0                    0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0              (.524+1j*.090)*.2];
        
load_power=1e3*[0           230+1j*142.5               0                230+1j*142.5          230+1j*142.5               0                     0                230+1j*142.5          230+1j*142.5               0                230+1j*142.5          137+1j*84              72+1j*45              72+1j*45              72+1j*45            13.5+1j*7.5            230+1j*142.5          230+1j*142.5          230+1j*142.5          230+1j*142.5          230+1j*142.5          230+1j*142.5          230+1j*142.5          230+1j*142.5          230+1j*142.5          230+1j*142.5          137+1j*85              75+1j*48              75+1j*48              75+1j*48              57+1j*34.5            57+1j*34.5            57+1j*34.5            57+1j*34.5];

source1=       [0                0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0        .1941*exp(1j*(-67.77)*pi/180)        0                     0                     0                     0                     0                     0                     0                     0                     0        .0702*exp(1j*(-124)*pi/180)        0                     0                     0                     0                     0                     0                     0         .1824*exp(1j*(-55.68)*pi/180)       0
                0                0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0        .1309*exp(1j*(11.9)*pi/180)          0                     0                     0                     0                     0                     0                     0                     0                     0        .0250*exp(1j*(-29.87)*pi/180)      0                     0                     0                     0                     0                     0                     0         .1190*exp(1j*(-84.11)*pi/180)       0
                0                0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0        .0758*exp(1j*(-7.13)*pi/180)         0                     0                     0                     0                     0                     0                     0                     0                     0        .0136*exp(1j*(-23.75)*pi/180)      0                     0                     0                     0                     0                     0                     0         .0573*exp(1j*(-143.56)*pi/180)      0
                0                0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0        .0586*exp(1j*(68.57)*pi/180)         0                     0                     0                     0                     0                     0                     0                     0                     0        .0075*exp(1j*(71.50)*pi/180)       0                     0                     0                     0                     0                     0                     0         .0401*exp(1j*(-175.58)*pi/180)      0
                0                0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0                     0        .0379*exp(1j*(46.53)*pi/180)         0                     0                     0                     0                     0                     0                     0                     0                     0        .0062*exp(1j*(77.21)*pi/180)       0                     0                     0                     0                     0                     0                     0         .0193*exp(1j*(111.39)*pi/180)       0  ]; 
            

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V_base=11e3;
S_base=10e6;
Z_base=V_base^2/S_base;
load_power=load_power/S_base;
line=line/Z_base;

Size_line=size(line);
Size_load_power=size(load_power);

[line_number,xxx]=size(find(line~= 0));
[xxx,load_power_number]=size(find(load_power~= 0));


%%%%   randomize loads with constant power factor with Laplacian random number generator  %%%%%%%

%     one=ones(Size_load_power(1,2),1);  
%     rand1=(60*one-15000*(0-(0.002/sqrt(2.))*sign(rand(Size_load_power(1,2),1)-0.5*one).*log(1*one-2*abs((rand(Size_load_power(1,2),1))-0.5*one))))/100;
%     
%     for kk=1:Size_load_power(1,2)
%         if rand1(kk,1)<0
%             rand1(kk,1)=(60-15000*(0-(0.002/sqrt(2.))*sign(rand(1,1)-0.5).*log(1-2*abs((rand(1,1))-0.5))))/100;
%         end
%     end
%     
%     load_power=load_power.*rand1';

    
    %%%%%%%%  for BIBC   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
BIBC=zeros(line_number,Size_load_power(1,2)-1);

for i=1:line_number
    for j=1:Size_load_power(1,2)-1
        if line(i,j)~= 0
            B=j;
            if i==1 && j==1
                BIBC(i,j)=1;
            else
                BIBC(:,j)=BIBC(:,i-1);
                BIBC(B,j)=1;
            end
        end
    end
end


for i=2:Size_load_power(1,2)
    
    if load_power(1,i)== 0
        BIBC(:,i-1)=0;
    end
end

%%%%%%%%      for BCBV      %%%%%%%%%%%%%%%%%%%%%%%%

BCBV=zeros(line_number,Size_load_power(1,2)-1);

            
for i=1:line_number
    for j=1:Size_load_power(1,2)-1
        if line(i,j)~= 0
            B=j;
            if i==1 && j==1
                BCBV(i,j)=line(i,j);
            else
                BCBV(j,:)=BCBV(i-1,:);
                BCBV(j,j)=line(i,j);
            end
        end
    end
end
 



% %%%%%%%%%%%%%   main  load flow        %%%%%%%%%%%%%%%% 
    
 DLF=BCBV*BIBC;
 V1=ones(Size_load_power(1,2)-1,1);
 V_bus=ones(Size_load_power(1,2)-1,1);
 I=zeros(Size_load_power(1,2)-1,1);
 Iaa=zeros(Size_load_power(1,2)-1,1);
 i=1;
 telorance=1;
 
 
 while i<=200   % maximum iterations
     for j=1:Size_load_power(1,2)-1
         I(j,1)=conj(load_power(1,j+1)/V_bus(j,1));
     end
     I_test(:,i)=I;
     V_bus= V1-(DLF*I);
     V_test(:,i)=V_bus;
     
     
     if i<=1
         telorance= 1;   % convergence condition
         
     else
         telorance= abs(abs(I_test(line_number,i))-abs(I_test(line_number,i-1)));
         
     
     end
     if abs(telorance) <= 1e-5
         fprintf('Power flow sloution found for %gth sample in "%g" iterations\n',jj,i)
         break
     end
        i=i+1;
 end
 if i==201
     fprintf('No soulotion, the algorithm is not converge\n')
     break
 end
 
 V_bus_size=size(V_bus);
 v_bus_shift=zeros(V_bus_size(1,1)+1,1);
 v_bus_shift(1,1)=1;
 
for i=2: V_bus_size(1,1)+1
    v_bus_shift(i,1)=V_bus(i-1,1);
end

 V_bus=v_bus_shift;
 
 I_bus_shift=zeros(V_bus_size(1,1)+1,1);
 I_bus_shift(1,1)=(V_bus(1,1)-V_bus(2,1))/line(1,1);
  
 for i=2: V_bus_size(1,1)+1
    I_bus_shift(i,1)=I(i-1,1);
 end

 I_bus=I_bus_shift;
 
%%%%%%  calculate the loads Impedance)

V_bus_size=size(V_bus);
Z_loads=zeros(Size_load_power(1,2),1);

for i=1:V_bus_size(1,1)
    if load_power(1,i) ~= 0
        Z_loads(i,1)=(V_bus(i,1))/conj(load_power(1,i)/V_bus(i,1));
    else
        Z_loads(i,1)=0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  harmonic load flow   %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=1;
line_1=real(line)+1j*imag(line)*h;
Z_loads_1=real(Z_loads)+1j*imag(Z_loads)*2*pi*50*h;
load_1=Z_loads_1';

for mm=1:5
    if mm==1
        h=5;
        source=abs(source1(mm,:)).*I_bus';
    end
    if mm==2
        h=7;
        source=abs(source1(mm,:)).*I_bus';
    end
    if mm==3
        h=11;
        source=abs(source1(mm,:)).*I_bus';
    end
    if mm==4
        h=13;
        source=abs(source1(mm,:)).*I_bus';
    end
    if mm==5
        h=17;
        source=abs(source1(mm,:)).*I_bus';
    end

 
    line=(real(line_1)+1j*imag(line_1)*h);
    Z_loads=(real(Z_loads_1)+1j*imag(Z_loads_1)*h);
    load=Z_loads';



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Size_line=size(line);
Size_load=size(load);
Size_source=size(source);

[line_number,xxx]=size(find(line~= 0));
[xxx,load_number]=size(find(load~= 0));
[xxx,source_number]=size(find(source~= 0));
parallel_number=source_number+load_number;

%%%%%%%%  for A   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A=zeros(line_number,parallel_number);
A1=A;
load_number1=load_number;
source_number1=source_number;
for i=Size_line:-1:1
    for j=Size_line:-1:1
        if line(i,j)~= 0
            if source(1,j+1)~= 0 
                A1(i,source_number1)=1;
                if source_number1 ~= 1 
                    source_number1=source_number1-1;
                end
            end
                if load(1,j+1)~= 0
                    A1(i,load_number1+source_number)=1;
                    load_number1=load_number1-1;
                end
        end
    end                    
end


for i=Size_line:-1:2
    for j=Size_line:-1:1
        if line(i,j)~= 0 
           for k=1:Size_line
               if line(k,i-1)~= 0               %% means these 2 line are conected to each other
                  for m=1:parallel_number
                      if A1(j,m)~= 0
                         A1(i-1,m)=A1(j,m);
                      end
                  end
               end
           end
        end
    end
end
 
A=A1;

%%%%%%%%%%%%%%%%%    for HA   %%%%%%%%%%%%%%%%%%%%%%

HA=zeros(line_number,parallel_number);
nn=0;
HA(1,:)=line(1,1);

for i=2:1:Size_line
    for j=1:1:Size_line
        if line(i,j)~= 0 
           for k=1:Size_line
               if line(k,i-1)~= 0
                  HA(j,:)=HA(i-1,:)+line(i,j)*A(j,:);
               end
           end
        end
    end
end

%%%%%%%%%%%%%%%%   for HAss   %%%%%%%%%%%%%%%%%%%
HAss=zeros(load_number);
i=1;
j=1;
for n=1:Size_load(1,2)
    if load(1,n)~= 0
        while j<=load_number
            HAss(i,j)=HA(n-1,source_number+j);
            j=j+1;
            
        end
        i=i+1;
        j=1;
    end
end
  
%%%%%%%%%%%%%  for HAsh    %%%%%%%%%%%%%%%%%%%%%

HAsh=zeros(load_number,source_number);
i=1;
j=1;

for n=1:Size_load(1,2)
    if load(1,n)~= 0
        while j<=source_number
            HAsh(i,j)=HA(n-1,j);
            j=j+1;
        end
        i=i+1;
        j=1;
    end
end
    
%%%%%%%%%%%%%%%    for Zs and HLF   %%%%%%%%%%%%%%%%%%%%%%%%%

Zs=zeros(load_number);
i=1;

for n=1:Size_load(1,2)
    if load(1,n)~= 0
        Zs(i,i)=load(1,n);
        i=i+1;
    end
end

HLF=HAss+Zs;

%%%%%%%%%%%%%%%%   finding currents of parallel elements  %%%%%%
Ih=zeros(source_number,1);
i=1;

for n=1:Size_source(1,2)
    if source(1,n)~= 0
        Ih(i,1)=source(1,n);
        i=i+1;
    end
end

Is=(HLF)^-1 * -HAsh*Ih;

%%%%%%%%%%%%  Harmonic load flow    %%%%%%%%%%%%%%%%%%%%%%%%%%

I=zeros(parallel_number,1);
size_Ih=size(Ih);
size_Is=size(Is);
k=1;
for i=1:parallel_number
    if i<= size_Ih(1,1)
        I(i,1)=Ih(i,1);
    else
        if k<= load_number
           I(i,1)=Is(k,1);
           k=k+1;
        end
    end
end

V_bus_h=HA*I;

 V_bus_h_size=size(V_bus_h);
 v_bus_h_shift=zeros(V_bus_h_size(1,1)+1,1);
 v_bus_h_shift(1,1)=0;
 
for i=2: V_bus_h_size(1,1)+1
    v_bus_h_shift(i,1)=V_bus_h(i-1,1);
end

 V_bus_h=v_bus_h_shift;

V_mag=abs(V_bus_h);
V_ang=angle(V_bus_h)*180/pi;
V_size=size(V_bus_h);
for i=1:V_size(1,1)
% fprintf('V%g=   %g , %g\n',i,V_mag(i,1),V_ang(i,1));
end

%%%%%%%%   injected current to each bus  %%%%%%%%%

for i=1:V_size(1,1)
    if source(1,i)~= 0
    I_harmonic_load(i,1)=V_bus_h(i,:)/Z_loads(i,1);
    else
    I_harmonic_load(i,1)=0;
    end
end

I_ingected=source'-I_harmonic_load; 

if mm==1
    V_bus_h_total_abs_5th(jj,:)=abs(V_bus_h');
    I_ingected_total_abs_5th(jj,:)=abs(I_ingected');
end
if mm==2
    V_bus_h_total_abs_7th(jj,:)=abs(V_bus_h');
    I_ingected_total_abs_7th(jj,:)=abs(I_ingected');
end
if mm==3
    V_bus_h_total_abs_11th(jj,:)=abs(V_bus_h');
    I_ingected_total_abs_11th(jj,:)=abs(I_ingected');
end
if mm==4
    V_bus_h_total_abs_13th(jj,:)=abs(V_bus_h');
    I_ingected_total_abs_13th(jj,:)=abs(I_ingected');
end
if mm==5
    V_bus_h_total_abs_17th(jj,:)=abs(V_bus_h');
    I_ingected_total_abs_17th(jj,:)=abs(I_ingected');
end

V_bus_h_total(jj+(mm-1)*sample_number,:)=abs(V_bus_h');
I_ingected_total(jj+(mm-1)*sample_number,:)=abs(I_ingected');

fprintf('sample No. %g for %gth harmonic generated.\n',jj,h);

end
 

end

fprintf('    ****   All  done   ****.\n');

V_bus_h_total_abs=abs(V_bus_h_total);
I_ingected_h_total_abs=abs(I_ingected_total);




