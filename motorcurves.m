function [maxPower, selTorque, limitation, gear] = motorcurves(jarda,vel)

if vel > 100
    vel = 100;
else
    %all is well
end

global final_ratio gbr motorrpm motortorque wheel_radius motors

motorvel = gbr.*(final_ratio*(vel/wheel_radius)); %engine rotational speed in rad/s
rpm = motorvel.*(60/(2*pi));                     %engine rpm
rpmmax = max(motorrpm);                          %max engine rpm
rpmmin = min(motorrpm);                          %min engine rpm

Torque = 0.98*0.98*jarda*interp1(motorrpm, motortorque, max(min(rpm,rpmmax),rpmmin));  %engine torque vector
Power = motors*Torque.*motorvel;                                                    %engine power vector
[maxPower, gear] = max(abs(Power));
selTorque = motors*Torque(gear);                                                    %selected engine torque

limitation = min(maxPower/vel);%, (selTorque*gbr(gear)*final_ratio)/wheel_radius);  %motor limitation