function [x_ecef,y_ecef,z_ecef,speed,xVel,yVel,zVel] = loadMammothData(path)
%% Import data from spreadsheet
% Script for importing data from the following spreadsheet:
%
%    Workbook: /Users/leila/Desktop/GPS Project/GNSS Raw Data Mammoth/Spreadsheets/Run1_03_01.xlsx
%    Worksheet: UsableData
%
% To extend the code for use with different selected data or a different
% spreadsheet, generate a function instead of a script.

% Auto-generated by MATLAB on 2020/03/21 15:54:54

%% Import the data
[~, ~, raw] = xlsread(path,'UsableData');
raw = raw(2:end,:);

%% Create output variable
data = reshape([raw{:}],size(raw));

%% Create table
Run10301 = table;

%% Allocate imported array to column variable names
x_ecef = data(:,1);
y_ecef = data(:,2);
z_ecef = data(:,3);
speed = data(:,4);
xVel = data(:,5);
yVel = data(:,6);
zVel = data(:,7);

%% Clear temporary variables
clearvars data raw;
end