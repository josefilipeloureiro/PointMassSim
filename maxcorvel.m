function corvel = TESTEmaxcorvelSS(raiocurva)

global by my
global mass weight
global SCZ ro

%Condições iniciais
Fz = weight;
lat_fric = my*(Fz/4) + by;
tiresy_lim = lat_fric*Fz;
corvel = sqrt((tiresy_lim/mass)*raiocurva);

velvector(1) = 0;
i = 2;
%300 iteracções
while i<300
    Fl = 0.5*ro*SCZ*corvel^2;
    Fz = (weight) + Fl;
    lat_fric = my*(Fz/4) + by; %muy = m*Fz + b;
    tiresy_lim = lat_fric*Fz;
    
    Fy = tiresy_lim;
    
    newvel = sqrt((Fy/mass)*raiocurva);
    velvector(i) = newvel;
    i = i+1;
    corvel = newvel;
end

%exponential moving average (útil para os casos em que oscila e não
%converge)
tenpct = round(length(velvector)*0.1);
velvector(1:tenpct)=[];
SMA = mean(velvector);
EMA(1) = SMA;
multiplier = 2/(length(velvector + 1));
for i = 1:length(velvector)
    EMA(i+1) = (velvector(i) - EMA(i))*multiplier + EMA(i);
end
corvel = EMA(end);