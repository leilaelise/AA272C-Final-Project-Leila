%% Load and Set Up Data

clear all;clc; close all;

% path = '/Users/leila/Desktop/GPS Project/Spreadsheets/Run1_2-01.xlsx'; 
% cutoff = 122;
%  path = '/Users/leila/Desktop/GPS Project/Spreadsheets/Run2_2-01.xlsx';
% cutoff = 1;
% path = '/Users/leila/Desktop/GPS Project/Spreadsheets/Run3_2-01.xlsx';
% cutoff = 1;
path = '/Users/leila/Desktop/GPS Project/Spreadsheets/Run4_2-01.xlsx'; %Downhill starts at x=27 
cutoff = 27;
% path = '/Users/leila/Desktop/GPS Project/Spreadsheets/Run5_2-01.xlsx';
% cutoff = 1; %Max Speed: 51.9431789320703
%  path = '/Users/leila/Desktop/GPS Project/Spreadsheets/Run6_2-01.xlsx';
%  cutoff = 1;
% path = '/Users/leila/Desktop/GPS Project/Spreadsheets/Run7_2-01.xlsx';
% cutoff = 120;
% path = '/Users/leila/Desktop/GPS Project/Spreadsheets/Run8_2-01.xlsx';
% cutoff = 198;
% path = '/Users/leila/Desktop/GPS Project/Spreadsheets/Run9_2-01.xlsx';
% cutoff = 1;
% path = '/Users/leila/Desktop/GPS Project/Spreadsheets/Run1_2-02.xlsx';
% cutoff = 92; %Max Speed: 39.2491528250011
% path = '/Users/leila/Desktop/GPS Project/Spreadsheets/Run2_2-02.xlsx';
% cutoff = 427;

[lat,long,speed,alt,time,svid] = loadRunData(path);

%% Get all data from GPS satellites

gpsLat = [];
gpsLong = [];
gpsSpeed = [];
gpsAlt = [];
gpsTime = [];


for i=1:length(svid)
    if (svid(i)=='GPS')
        gpsLat = [gpsLat lat(i)];
        gpsLong = [gpsLong long(i)];
        gpsSpeed = [gpsSpeed speed(i)];
        gpsAlt = [gpsAlt alt(i)];
        gpsTime = [gpsTime time(i)];
    end
end

gpsLat = gpsLat';
gpsLong = gpsLong';
gpsSpeed = gpsSpeed';
gpsAlt = gpsAlt';
gpsTime = gpsTime'/10^3; %Convert to [s]

gpsData = [gpsLat gpsLong gpsSpeed gpsAlt gpsTime];

%% Set up measurement vector

ecef = lla2ecef([gpsLat gpsLong gpsAlt]);
x_meas = ecef(:,1);
y_meas = ecef(:,2);
z_meas = ecef(:,3);

thetaMeas = atan2(y_meas,x_meas);
phiMeas = atan2(z_meas,sqrt(x_meas.^2+y_meas.^2));
%zGPS = [x y z theta phi v];
meas = [x_meas y_meas z_meas thetaMeas phiMeas gpsSpeed];
mu = [x_meas(1) y_meas(1) z_meas(1) thetaMeas(1) phiMeas(1) gpsSpeed(1)]';
%% Set Up Kalman Filter

n = 6;
m = 6;

dt = 1;

H = eye(6);
P = eye(6);
%TODO: use actual computed standard deviations from phone and tune R and Q
sigmaMeas = 10; %[m]
sigmaProcess = 0.05*10^3; %[m]
R = eye(6)*sigmaMeas^2;
Q = eye(6)*sigmaProcess^2;

latTraj = [];
longTraj = [];
altTraj = [];
speedTraj = [];

for t = 2:length(gpsTime)
    
    %Predict
    F = [1 0 0 0 0 cos(phiMeas(t))*cos(thetaMeas(t))*dt; 0 1 0 0 0 cos(phiMeas(t))*sin(thetaMeas(t))*dt; 0 0 1 0 0 sin(phiMeas(t))*dt;...
        0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1];
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
    
    P = (eye(6) - Kt*H)*P*transpose(eye(6)-Kt*H)+Kt*R*transpose(Kt);
    yt = meas(t,:)-H*mu;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the trajectory of the agent
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

latTraj = latTraj(cutoff:end);
longTraj = longTraj(cutoff:end);
altTraj = altTraj(cutoff:end);
speedTraj = speedTraj(cutoff:end);

outputVec = [latTraj' longTraj' speedTraj'];

% figure
plot3(longTraj,latTraj,altTraj)
grid on
xlabel('Longitude');
ylabel('Latitude');
zlabel('Altitude');


% hold on
% p = plot3(longTraj(1),latTraj(1),altTraj(1),'o','MarkerFaceColor','red');
% hold off
% axis manual
% 
% for k = 2:length(longTraj)
%     p.XData = longTraj(k);
%     p.YData = latTraj(k);
%     p.ZData = altTraj(k);
%     drawnow, pause(0.01*gpsSpeed(k))
% end
% 
%  axis([min(longTraj) max(longTraj) min(latTraj) max(latTraj) min(altTraj) max(altTraj)])
% 
figure
subplot(2,1,1)
plot(longTraj,latTraj)
xlabel('Longitude');
ylabel('Latitude');
grid on
% hold on
% p = plot(longTraj(1),latTraj(1),'o','MarkerFaceColor','red');
% hold off
% axis manual
% 
% for k = 2:length(longTraj)
%     p.XData = longTraj(k);
%     p.YData = latTraj(k);
%     drawnow
%     pause(0.01*gpsSpeed(k))
% end
subplot(2,1,2)
plot(1:length(altTraj), altTraj)
xlabel('Time[s]');
ylabel('Altitude[m]');
grid on

max(speedTraj)*2.23694

%% Save data for Google Maps API
dlmwrite('lat.csv', latTraj, 'delimiter', ',', 'precision', 12);
dlmwrite('long.csv', longTraj, 'delimiter', ',', 'precision', 12);
dlmwrite('speed.csv', speedTraj, 'delimiter', ',', 'precision', 12);
