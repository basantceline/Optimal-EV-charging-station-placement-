clc;
clear all;

m = load('loaddata.m');
l = load('linedata.m');
voltageBase = indLoadFlow(m, l);
total = sum(m(:,2));
for i= 2:28
% Choose bus to analyze
bus_to_analyze = i;

% Initialize arrays
x = [];
voltage12 = [];

% Set maximum load increment (e.g., 300% of base load)
max_load_increment = 3 * m(bus_to_analyze, 2);

% Set smaller step size
step_size = max_load_increment / 100;

m1 = m;
for increment = 0:step_size:max_load_increment
    m1(bus_to_analyze, 2) = m(bus_to_analyze, 2) + increment;
    
    voltagePerturbed = indLoadFlow(m1, l);
    
    if isempty(voltagePerturbed)
        disp(['Load flow did not converge at increment: ', num2str(increment)]);
        break;
    end
    
    x(end+1) = m1(bus_to_analyze, 2) / m(bus_to_analyze, 2);  % Load as multiple of base load
    voltage12(end+1) = voltagePerturbed(bus_to_analyze);
end

if i==2 
   figure;
plot(x, voltage12, 'b-', 'LineWidth', 2);
xlabel('Load (multiple of base load)');
ylabel('Voltage (p.u.)');
title(['PV Curve for Bus ', num2str(bus_to_analyze)]);
grid on;
xlim([1 max(x)]);
ylim([min(voltage12)-0.05 max(voltage12)+0.05]);
else
    plot(x, voltage12, 'b-', 'LineWidth', 2);
end

    
end
    %      delV=abs((voltageBase(i)-voltagePerturbed(i))*11000)
    %      delP=abs(m(i,2)-m1(i,2))
    %      slope(1,end+1)=delV/delP;

% dataActual=table();
% for i= 2:29
%     m=load('loaddata.m');
%     l=load('linedata.m');
%     data=m(i,2)+120;
%     m(i,2)= data;
%     voltage= indLoadFlow(m,l);
%     if i == 2
%         dataActual.Pentration2= voltage
%     else
%          name = "Pentration at "+ num2str(i)
%         dataActual.(name)= voltage;
%     end
%
% end
% names= strings(1, 29);
% for q= 1:29
% names(1,q)="Bus "+num2str(q);
% end
% dataActual.Properties.RowNames=names;
% writetable(dataActual,'voltageProfile.xlsx','WriteRowNames',true)




function voltageModify = indLoadFlow(m,l)
format short;
br=length(l);
no=length(m);
f=0;
d=0;
MVAb=100;
KVb=11;
Zb=(KVb^2)/MVAb;
Pg = zeros(no,1);
Pg1 = zeros(no,1);

for i=1:br
    R(i,1)=(l(i,4))/Zb;
    X(i,1)=(l(i,5))/Zb;
end
for i=1:no
    P(i,1)=(m(i,2)/(1000*MVAb));
    Q(i,1)=(m(i,3)/(1000*MVAb));
end
R;
X;
P;
Q;

C=zeros(br,no);
for i=1:br
    a=l(i,2);
    b=l(i,3);
    for j=1:no
        if a==j
            C(i,j)=-1;
        end
        if b==j
            C(i,j)=1;
        end
    end
end
C;
e=1;
for i=1:no
    d=0;
    for j=1:br
        if C(j,i)==-1
            d=1;
        end
    end
    if d==0
        endnode(e,1)=i;
        e=e+1;
    end
end
endnode;
h=length(endnode);
for j=1:h
    e=2;

    f=endnode(j,1);
    % while (f~=1)
    for s=1:no
        if (f~=1)
            k=1;
            for i=1:br
                if ((C(i,f)==1)&&(k==1))
                    f=i;
                    k=2;
                end
            end
            k=1;
            for i=1:no
                if ((C(f,i)==-1)&&(k==1));
                    f=i;
                    g(j,e)=i;
                    e=e+1;
                    k=3;
                end
            end
        end
    end
end
for i=1:h
    g(i,1)=endnode(i,1);
end
g;
w=length(g(1,:));
for i=1:h
    j=1;
    for k=1:no
        for t=1:w
            if g(i,t)==k
                g(i,t)=g(i,j);
                g(i,j)=k;
                j=j+1;
            end
        end
    end
end
g;
for k=1:br
    e=1;
    for i=1:h
        for j=1:w-1
            if (g(i,j)==k)
                if g(i,j+1)~=0
                    adjb(k,e)=g(i,j+1);
                    e=e+1;
                else
                    adjb(k,1)=0;
                end
            end
        end
    end
end
adjb;
for i=1:br-1
    for j=h:-1:1
        for k=j:-1:2
            if adjb(i,j)==adjb(i,k-1)
                adjb(i,j)=0;
            end
        end
    end
end
adjb;
x=length(adjb(:,1));
ab=length(adjb(1,:));
for i=1:x
    for j=1:ab
        if adjb(i,j)==0 && j~=ab
            if adjb(i,j+1)~=0
                adjb(i,j)=adjb(i,j+1);
                adjb(i,j+1)=0;
            end
        end
        if adjb(i,j)~=0
            adjb(i,j)=adjb(i,j)-1;
        end
    end
end
adjb;
for i=1:x-1
    for j=1:ab
        adjcb(i,j)=adjb(i+1,j);
    end
end
b=length(adjcb);

% voltage current program

for i=1:no
    vb(i,1)=1;
end
for s=1:10
    for i=1:no
        nlc(i,1)=conj(complex(P(i,1),Q(i,1)))/(vb(i,1));
    end
    nlc;
    for i=1:br
        Ibr(i,1)=nlc(i+1,1);
    end
    Ibr;
    xy=length(adjcb(1,:));
    for i=br-1:-1:1
        for k=1:xy
            if adjcb(i,k)~=0
                u=adjcb(i,k);
                Ibr(i,1)=Ibr(i,1)+Ibr(u,1);
            end
        end
    end
    Ibr;
    for i=2:no
        g=0;
        for a=1:b
            if xy>1
                if adjcb(a,2)==i-1
                    u=adjcb(a,1);
                    vb(i,1)=((vb(u,1))-((Ibr(i-1,1))*(complex((R(i-1,1)),X(i-1,1)))));
                    g=1;
                end
                if adjcb(a,3)==i-1
                    u=adjcb(a,1);
                    vb(i,1)=((vb(u,1))-((Ibr(i-1,1))*(complex((R(i-1,1)),X(i-1,1)))));
                    g=1;
                end
            end
        end
        if g==0
            vb(i,1)=((vb(i-1,1))-((Ibr(i-1,1))*(complex((R(i-1,1)),X(i-1,1)))));
        end
    end
    s=s+1;
end
nlc;
Ibr;
vb;
vbp=[abs(vb)];

for i=1:no
    va(i,2)=vbp(i,1);
end
for i=1:no
    va(i,1)=i;
    P1(i) = P(i);
    Q1(i) = Q(i);
end

va;
Ibrp=[abs(Ibr)];
PL(1,1)=0;
QL(1,1)=0;

% losses at base case
for f=1:br
    Pl(f,1)=(Ibrp(f,1)^2)*R(f,1);
    Ql(f,1)=X(f,1)*(Ibrp(f,1)^2);
    PL(1,1)=PL(1,1)+Pl(f,1);
    QL(1,1)=QL(1,1)+Ql(f,1);
end

Plosskw=(Pl)*100000;
Qlosskw=(Ql)*100000;
PL=(PL)*100000;
QL=(QL)*100000;

voltage = vbp(:,1);

v_mag = va(:,2);
voltageModify=v_mag;
% % for plotting bar and formatting the graph
% bus=1:1:29;
% bus=bus.';
% bar(bus,voltage,0.2)
% xticks(bus);
% xlabel('Bus');
% ylabel('Voltage in pu');
% ylim([0 1.1]);
% % for plotting bar and formatting the graph
% bus=1:1:no;
% bus=bus.';
% pBus= bus(2:end);
% %subplot(2,1,1)
% bar(pBus,Plosskw);
% xticks(pBus);
% xticklabels({'1-2','2-3','3-4','4-5','5-6','6-7','7-8','8-9','9-10','10-11','11-12','12-13','13-14','2-15','5-16','16-17','8-18','18-19','19-20','20-21','11-22','22-23','22-24','24-25','12-26','26-27','26-28','28-29'})
% xlabel('Branch');
% ylabel('Active Power loss (Kw)');
% %ylim([0 1.1]);
% %subplot(2,1,2)
% figure()
% bar(bus,voltage,0.5);
% xticks(bus);
% xlabel('Bus');
% ylabel('Voltage in pu');
% ylim([0 1.1]);
% figure()
% bar(pBus,Qlosskw);
% xticks(pBus);
% xticklabels({'1-2','2-3','3-4','4-5','5-6','6-7','7-8','8-9','9-10','10-11','11-12','12-13','13-14','2-15','5-16','16-17','8-18','18-19','19-20','20-21','11-22','22-23','22-24','24-25','12-26','26-27','26-28','28-29'})
% xlabel('Branch');
% ylabel('Reactive Power loss (KVAR)');
end


