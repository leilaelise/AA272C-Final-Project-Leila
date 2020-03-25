%ProcessGnssMeasScript.m, script to read GnssLogger output, compute and plot:
clear all; clc;
% pseudoranges, C/No, and weighted least squares PVT solution
%
% you can run the data in pseudoranges log files provided for you:
prFileName = 'gnss_log_2020_03_01_10_30_47.txt'; %with duty cycling, no carrier phase
fileName = 'Run2_3_01.xlsx';
% prFileName = 'pseudoranges_log_2016_08_22_14_45_50.txt'; %no duty cycling, with carrier phase
% as follows
% 1) copy everything from GitHub google/gps-measurement-tools/ to
%    a local directory on your machine
% 2) change 'dirName = ...' to match the local directory you are using:
dirName = '~/Desktop/GPS Project/GNSS Processing/gps tools/opensource/Ski Files/Raw Data';
% 3) run ProcessGnssMeasScript.m script file
param.llaTrueDegDegM = [];

%Author: Frank van Diggelen
%Open Source code for processing Android GNSS Measurements

%% data
%To add your own data:
% save data from GnssLogger App, and edit dirName and prFileName appropriately
%dirName = 'put the full path for your directory here';
%prFileName = 'put the pseuoranges log file name here';

%% parameters
%param.llaTrueDegDegM = [];
%enter true WGS84 lla, if you know it:
param.llaTrueDegDegM = [37.422578, -122.081678, -28];%Charleston Park Test Site

%% Set the data filter and Read log file
dataFilter = SetDataFilter;
[gnssRaw,gnssAnalysis] = ReadGnssLogger(dirName,prFileName,dataFilter);
if isempty(gnssRaw), return, end

%% Get online ephemeris from Nasa ftp, first compute UTC Time from gnssRaw:
fctSeconds = 1e-3*double(gnssRaw.allRxMillis(end));
utcTime = Gps2Utc([],fctSeconds);
allGpsEph = GetNasaHourlyEphemeris(utcTime,dirName);
if isempty(allGpsEph), return, end

%% process raw measurements, compute pseudoranges:
[gnssMeas] = ProcessGnssMeas(gnssRaw);

%% plot pseudoranges and pseudorange rates
h1 = figure;
[colors] = PlotPseudoranges(gnssMeas,prFileName);
%%

h2 = figure;
PlotPseudorangeRates(gnssMeas,prFileName,colors);
h3 = figure;
PlotCno(gnssMeas,prFileName,colors);

%% compute WLS position and velocity
gpsPvt = GpsWlsPvt(gnssMeas,allGpsEph);

%% plot Pvt results
h4 = figure;
ts = 'Raw Pseudoranges, Weighted Least Squares solution';
PlotPvt(gpsPvt,prFileName,param.llaTrueDegDegM,ts); drawnow;
h5 = figure;
PlotPvtStates(gpsPvt,prFileName);

%% Plot Accumulated Delta Range
if any(any(isfinite(gnssMeas.AdrM) & gnssMeas.AdrM~=0))
    [gnssMeas]= ProcessAdr(gnssMeas);
    h6 = figure;
    PlotAdr(gnssMeas,prFileName,colors);
    [adrResid]= GpsAdrResiduals(gnssMeas,allGpsEph,param.llaTrueDegDegM);drawnow
    h7 = figure;
    PlotAdrResids(adrResid,gnssMeas,prFileName,colors);
end
%% Convert results to x,y,z ECEF coordinates
%Remove NaN values
pos = gpsPvt.allLlaDegDegM;
vel = gpsPvt.allVelMps;
A = pos(:,1);
B = pos(:,2);
C = pos(:,3);
lat = A(~isnan(A));
long = B(~isnan(B));
alt = C(~isnan(B));
posLLA = [lat long alt];
A = vel(:,1);
B = vel(:,2);
C = vel(:,3);
xVel = A(~isnan(A));
yVel = B(~isnan(B));
zVel = C(~isnan(B));

%Convert to ECEF
xECEF = [];
yECEF = [];
zECEF = [];
xDotECEF = [];
yDotECEF = [];
zDotECEF = [];
for i=1:length(posLLA)
    posECEF = lla2ecef(posLLA(i,:));
    xECEF = [xECEF posECEF(1)];
    yECEF = [yECEF posECEF(2)];
    zECEF = [zECEF posECEF(3)];
    [xDot,yDot,zDot] = ned2ecefv(xVel(i),yVel(i),zVel(i),lat(i),long(i));
    xDotECEF = [xDotECEF xDot];
    yDotECEF = [yDotECEF yDot];
    zDotECEF = [zDotECEF zDot];
end

xECEF = xECEF';
yECEF = yECEF';
zECEF = zECEF';
xDotECEF = xDotECEF';
yDotECEF = yDotECEF';
zDotECEF = zDotECEF';

vel = [];
for i=1:length(xDotECEF)
vel = [vel norm([xDotECEF(i) yDotECEF(i) zDotECEF(i)])];
end

state = [xECEF yECEF zECEF xDotECEF yDotECEF zDotECEF];
speed = vel';
stateV2 = [xECEF yECEF zECEF speed];
time = [0:length(state)-1]';


T = table(xECEF, yECEF, zECEF, speed, xDotECEF, yDotECEF, zDotECEF);
cd '~/Desktop/GPS Project/GNSS Raw Data Mammoth/Spreadsheets';
writetable(T,fileName,'Sheet','UsableData','Range','A1')

%% end of ProcessGsnssMeasScript
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2016 Google Inc.
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%     http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
