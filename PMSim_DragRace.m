%% Load Track
[TrackFileName,TrackFilePath] = uigetfile('*.cir','Select the track');
fid = fopen(TrackFileName);
i = 1;
while ~feof(fid)
    linha = fgetl(fid);
    linha = str2num(linha);
    pos(i) = linha(1);
    raio_local(i) = linha(2);
    phi(i) = linha(3);
    i = i+1;
end
fclose(fid);
clear i
%% Init
global final_ratio gbr motorrpm motortorque motors
global by my bx mx
global mass weight WD_front wheel_radius
global SCX SCZ ro

%Load Car
[VecFileName,VecFilePath] = uigetfile('*_vec.mat','Select the vehicle');
load(VecFileName)
ro = 101325/(287.058*293.15); %Densidade do ar a 20ºC e 1atm

%Load Tyres
[TyresFileName,TyresFilePath] = uigetfile('*_tyres.mat','Select the tyres');
load(TyresFileName)

%Load Track
%load('') load já foi na secção anterior

%Number of segments
seg = length(pos);

raio_local = abs(raio_local);   %radius of each segment
d = zeros(seg,1);               %distance of each segment

for i = 1:length(pos)
    if raio_local(i) > 0.10e+04
        raio_local(i) = 0.10e+04;
    else 
        %all is well
    end
    
    if i < length(pos)
        d(i) = pos(i+1) - pos(i);
    elseif i == length(pos)
        d(i) = pos(i) - pos(i-1);
    else
        fprintf('Erro no cálculo das distâncias de segmentos.')
    end
end

jarda = 1;%input('Torque factor, a.k.a. jarda: '); %[0;1]

checkeredVel = 0;

%Initial Conditions
vel = zeros(seg,1); vel(1) = checkeredVel;
Fd = zeros(seg,1); Fd(1) = 0.5*ro*SCX*vel(1)^2;
Fl = zeros(seg,1); Fl(1) = 0.5*ro*SCZ*vel(1)^2;
Fz = zeros(seg,1); Fz(1) = weight + Fl(1);
motor_lim = zeros(seg,1); [~, ~, motor_lim(1), ~] = motorcurves(jarda,vel(1));
Fe = zeros(seg,1); Fe(1) = motor_lim(1);
accx = zeros(seg,1); accx(1) = (min((mx*(Fz(1)/4) + bx)*Fz(1)*(1-WD_front), motor_lim(1))-Fd(1)-Fz(1)*0.01)/mass;
accy = zeros(seg,1); accy(1) = 0;
time = zeros(seg,1); time(1) = (-vel(1)+sqrt((vel(1)^2)+2*accx(1)))/(accx(1));
elapsed_time = zeros(seg,1); elapsed_time(1) = time(1);

clear i
% Forward simulation (sem travagens)
for i = 2:length(pos)
    corvel = maxcorvel(raio_local(i)); %maximum corner velocity
    
    vel(i) = min(corvel,vel(i-1) + accx(i-1)*time(i-1));
    Fd(i) = 0.5*ro*SCX*vel(i)^2;
    Fl(i) = 0.5*ro*SCZ*vel(i)^2;
    Fz(i) = weight + Fl(i);
    lat_fric = my*(Fz(i)/4) + by; %muy = m*Fz + b;
    tiresy_lim = lat_fric*Fz(i);
    long_fric = mx*(Fz(i)*(1-WD_front)/4) + bx; %mux = m*Fz + b;
    tiresx_lim = long_fric*Fz(i)*(1-WD_front);

    [~, ~, limitation, ~] = motorcurves(jarda,vel(i));
    motor_lim = limitation;

    Fy = min(tiresy_lim, mass*(vel(i)^2)/raio_local(i));
    Fe(i) = min(sqrt(1 - (Fy/tiresy_lim)^2)*tiresx_lim, motor_lim);
    accx(i) = (Fe(i)-Fd(i)-Fz(i)*0.01)/mass;
    accy(i) = Fy/mass;

    time(i) = d(i)/vel(i);
    elapsed_time(i) = sum(time) + time(i);
end

clear i j k
%% Reverse Simulation (com travagens)
revvel = zeros(seg,1); revvel(1) = 80.4444;
revFd = zeros(seg,1); revFd(1) = 0.5*ro*SCX*revvel(1)^2;
revFl = zeros(seg,1); revFl(1) = 0.5*ro*SCZ*revvel(1)^2;
revFz = zeros(seg,1); revFz(1) = weight + revFl(1);
revaccx = zeros(seg,1); revaccx(1) = (((mx*(revFz(1)/4) + bx)*revFz(1))+revFd(1)+revFz(1)*0.01)/mass;
revaccy = zeros(seg,1); revaccy(1) = 0;
revtime = zeros(seg,1); revtime(1) = (-revvel(1)+sqrt((revvel(1)^2)+2*revaccx(1)))/(revaccx(1));
revelapsed_time = zeros(seg,1); revelapsed_time(1) = revtime(1);

j = 2;
for i = length(pos):-1:2
    revcorvel = maxcorvel(raio_local(i)); %maximum corner velocity

    revvel(j) = min(revcorvel, revvel(j-1) + revaccx(j-1)*revtime(j-1));
    revFd(j) = 0.5*ro*SCX*revvel(j)^2;
    revFl(j) = 0.5*ro*SCZ*revvel(j)^2;
    revFz(j) = weight + revFl(j);
    lat_fric = my*(Fz(j)/4) + by; %muy = m*Fz + b;
    tiresy_lim = lat_fric*Fz(j);
    long_fric = mx*(Fz(j)/4) + bx; %mux = m*Fz + b;
    tiresx_lim = long_fric*Fz(j);

    Fy = min(tiresy_lim, mass*(revvel(j)^2)/raio_local(i));
    revFx = sqrt(1 - (Fy/tiresy_lim)^2)*tiresx_lim;
    revaccx(j) = (revFx+revFd(j)+revFz(j)*0.01)/mass;
    revaccy(j) = Fy/mass;

    revtime(j) = d(j)/revvel(j);
    revelapsed_time(j) = sum(revtime) + revtime(j);
    
    j = j+1;
end

clear i j k
j = 1;
for i = length(revvel):-1:1
    newrevvel(j) = revvel(i);
    newrevFd(j) = revFd(i);
    newrevFz(j) = revFz(i);
    newrevaccx(j) = revaccx(i);
    newrevaccy(j) = revaccy(i);
    j = j+1;
end
clear i j k revvel revFd revFz revaccx revaccy
revvel = newrevvel;
revFd = newrevFd;
revFz = newrevFz;
revaccx = newrevaccx;
revaccy = newrevaccy;
%% Minimum of vel and revvel
Vel_PM = zeros(seg,1);
Time_PM = zeros(seg,1);
Deltatime_PM = zeros(seg,1);
Fd_PM = zeros(seg,1);
Fz_PM = zeros(seg,1);
Fe_PM = zeros(seg,1);
AccX_PM = zeros(seg,1);
AccY_PM = zeros(seg,1);
for i = 1:length(vel)
    Vel_PM(i) = min(vel(i), revvel(i));
    
    if Vel_PM(i) == vel(i)
        Fd_PM(i) = Fd(i);
        Fz_PM(i) = Fz(i);
        Fe_PM(i) = Fe(i);
        AccX_PM(i) = accx(i);
        AccY_PM(i) = accy(i);
    elseif Vel_PM(i) == revvel(i)
        Fd_PM(i) = revFd(i);
        Fz_PM(i) = revFz(i);
        AccX_PM(i) = -revaccx(i);
        AccY_PM(i) = revaccy(i);
    else
        fprintf('Erro a definir o perfil de carga vertical Fz ou aceleração longitudinal Accx.')
    end
    
    if Vel_PM(i) == 0
        Time_PM(i) = sqrt(2/accx(i));
    else
        Time_PM(i) = d(i)/Vel_PM(i);
    end
    Deltatime_PM(i) = sum(Time_PM) + Time_PM(i);
end
% figure
% plot(vel*3.6)
% hold on
% plot(revvel*3.6)
% hold off
% xlabel('Track Length (m)'), ylabel('Velocity (km/h)'), legend('Forward', 'Reverse')
% figure
% load('JKVEstoril.1.38.295.mat')
% plot(x/10, Velkmh)
% hold on
% plot(Vel_PM*3.6)
% xlabel('Track Length (m)'), ylabel('Velocity (km/h)')
% legend('Real', 'PMSim')
clear i j k
%% Energy Consumption
EnergyRR_PM = zeros(seg,1);
EnergyHC_PM = zeros(seg,1);
EnergyAero_PM = zeros(seg,1);
EnergyLI_PM = zeros(seg,1);
MechanicalEnergy = zeros(seg,1);
MechanicalEnergy_PM = zeros(seg,1);
for i = 1:seg
    %if Vel_PM(i) == vel(i)
    if AccX_PM(i) >= 0
        %MechanicalEnergy [kWh] = 2.778*10^-7 * [m*g(Crr*cos(phi) + sin(phi)) + 0.5*ro*Cx*A*v^2 + (m + (I*gr/eta*r_motor^2))*a] * d[m]
        EnergyRR_PM(i) = 2.778*10^(-7) * (Fz_PM(i)*0.01) * cos(phi(i)) * d(i);
        EnergyHC_PM(i) = 0;
        EnergyAero_PM(i) = 2.778*10^(-7) * 0.5 * (Vel_PM(i)^2) * ro * SCX * d(i);
        %EnergyLI_PM(i) = 2.778*10^(-7) * (Fz_PM(i)/9.81) * AccX_PM(i) * d(i);
        EnergyLI_PM(i) = 2.778*10^(-7) * Fe_PM(i) * d(i);
        %EnergyAI_PM(i) = (I*gr^2 / eta*r^2) * AccX_OK(i) * d(i);
        MechanicalEnergy(i) = EnergyRR_PM(i) + EnergyHC_PM(i) + EnergyLI_PM(i) + EnergyAero_PM(i);
        MechanicalEnergy_PM(i) = sum(MechanicalEnergy);
    end
end
%fprintf('Lap time: %f s\nEnergia gasta: %f kWh\nTop Speed: %f km/h\n', Deltatime_PM(end), MechanicalEnergy_PM(end), max(Vel_PM)*3.6)
Results = zeros(3,1);
Results(1) = Deltatime_PM(end);
Results(2) = MechanicalEnergy_PM(end);
Results(3) = max(Vel_PM)*3.6;
PMSimGUI(Vel_PM, Deltatime_PM, MechanicalEnergy_PM)
clearvars -except Vel_PM Time_PM Deltatime_PM Fd_PM Fz_PM Fe_PM AccX_PM AccY_PM EnergyRR_PM EnergyHC_PM EnergyAero_PM EnergyLI_PM MechanicalEnergy MechanicalEnergy_PM