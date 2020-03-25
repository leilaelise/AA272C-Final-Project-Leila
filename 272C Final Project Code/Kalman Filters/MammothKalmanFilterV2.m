%% Kalman Filter Using GNSS WLS solution data, Mammoth data

clear all;clc; close all;

%Data from the run that I fell
path = '/Users/leila/Desktop/GPS Project/GNSS Raw Data Mammoth/Spreadsheets/Run2_3_01.xlsx';
cutoff = 300;
endCutoff = 2048;

%Data from after fall, when Zach took the phone
% path = '/Users/leila/Desktop/GPS Project/GNSS Raw Data Mammoth/Spreadsheets/Run2_3_01.xlsx';
%cutoff = 1;
%endCutoff = 0;

%Data from before fall
%path = '/Users/leila/Desktop/GPS Project/GNSS Raw Data Mammoth/Spreadsheets/Run3_3_01.xlsx';
% cutoff = 300;
% endCutoff = 2048;

% path = '/Users/leila/Desktop/GPS Project/GNSS Raw Data Mammoth/Spreadsheets/Run4_3_01.xlsx';
% cutoff = 137;
% endCutoff = 0;
% 
% path = '/Users/leila/Desktop/GPS Project/GNSS Raw Data Mammoth/Spreadsheets/Zach1_3_01.xlsx';
% cutoff = 1;
% endCutoff = 0;



[x_meas,y_meas,z_meas,speed_meas] = loadMammothData(path);

%Remove anomalies
x_ecef = x_meas;
y_ecef = y_meas;
z_ecef = z_meas;
speed = speed_meas;
for i=1:length(speed_meas)
    if (abs(speed_meas(i))>30) %Max speed at 70 mph
        x_ecef(i) = x_ecef(i-1);
        y_ecef(i) = y_ecef(i-1);
        z_ecef(i) = z_ecef(i-1);
        speed(i) = speed(i-1);
    end
end


time = [1:length(speed)]';
%Set up initial matrices and dimensions
n = 6;

thetaMeas = atan2(y_ecef,x_ecef);
phiMeas = atan2(z_ecef,sqrt(x_ecef.^2+y_ecef.^2));
accel = [0 diff(speed)']';
meas = [x_ecef y_ecef z_ecef thetaMeas phiMeas speed accel];
mu = [x_ecef(1) y_ecef(1) z_ecef(1) thetaMeas(1) phiMeas(1) speed(1) accel(1)]';

%Set Up Kalman Filter
n = 7;

dt = 1;

H = eye(n);
P = eye(n);
%TODO: use actual computed standard deviations from phone and tune R and Q
stdX = 52.846; %[m]
stdY = 46.117; %[m]
stdZ = 65.507; %[m]
stdXDot = 1.392; %[m/s]
stdYDot = 1.183; %[m/s]
stdZDot = 1.514; %[m/s]
stdV = 2.37; %[m/s]
stdA = 2.44; %[m/s^2]

% R = [stdX^2 0 0 0 0 stdV^2 0; 0 stdY^2 0 0 0 stdV^2 0; 0 0 stdZ^2 0 0 stdV^2 0;...
%      0 0 0 stdX*stdY 0 stdV^2 0; 0 0 0 0 stdZ^2 stdV^2 0; 0 0 0 0 0 stdV^2 0; 0 0 0 0 0 stdV^2 stdA^2];
%R = eye(n);
load('Rcov.mat');

Q = eye(n)*0.15;

latTraj = [];
longTraj = [];
altTraj = [];
speedTraj = [];

for t = 2:length(time)
    
    %Predict
    F = [1 0 0 0 0 cos(phiMeas(t))*cos(thetaMeas(t))*dt 0.5*cos(phiMeas(t))*cos(thetaMeas(t))*dt^2;...
         0 1 0 0 0 cos(phiMeas(t))*sin(thetaMeas(t))*dt 0.5*cos(phiMeas(t))*sin(thetaMeas(t))*dt^2;...
         0 0 1 0 0 sin(phiMeas(t))*dt 0.5*sin(phiMeas(t))*dt^2;...
         0 0 0 1 0 0 dt; 0 0 0 0 1 0 dt; 0 0 0 0 0 1 dt; 0 0 0 0 0 0 1];
    mu = F*mu;
    P = F*P*transpose(F) + Q;
    
    %Update
    yt = meas(t,:)'-H*mu;
    Kt = P*transpose(H)*inv(R + H*P*transpose(H));
    mu = mu + Kt*yt;
    
    lla = ecef2lla([mu(1) mu(2) mu(3)]);
    latTraj = [latTraj lla(1)];
    longTraj = [longTraj lla(2)];
    altTraj = [altTraj lla(3)];
    speedTraj = [speedTraj mu(6)];
    
    P = (eye(n) - Kt*H)*P*transpose(eye(n)-Kt*H)+Kt*R*transpose(Kt);
    yt = meas(t,:)-H*mu;
end


if (endCutoff~=0)
    latTraj = latTraj(cutoff:endCutoff);
    longTraj = longTraj(cutoff:endCutoff);
    altTraj = altTraj(cutoff:endCutoff);
    speedTraj = speedTraj(cutoff:endCutoff);
else
    %No cutoff
    latTraj = latTraj(cutoff:end);
    longTraj = longTraj(cutoff:end);
    altTraj = altTraj(cutoff:end);
    speedTraj = speedTraj(cutoff:end);
end


traj = [latTraj' longTraj'];


%% 

figure
subplot(2,1,1)
plot(longTraj,latTraj)
xlabel('Longitude');
ylabel('Latitude');
grid on

subplot(2,1,2)
plot(1:length(altTraj), altTraj)
xlabel('Time[s]');
ylabel('Altitude[m]');
grid on
%% Save data for Google Maps API
dlmwrite('lat.csv', latTraj, 'delimiter', ',', 'precision', 12);
dlmwrite('long.csv', longTraj, 'delimiter', ',', 'precision', 12);
dlmwrite('speed.csv', speedTraj, 'delimiter', ',', 'precision', 12);