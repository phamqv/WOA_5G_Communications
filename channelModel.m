%==========================================================================
% Authors: Quoc-Viet PHAM
% Created: 11/10/2017
% Current: 04/05/2019
% E-mail:  vietpq90@gmail.com
%
% This function returns the channel gain matrix in linear, the row is the
% index of mobile users and the column is the index of SeNBs
%==========================================================================
function [hArray, dArray, eNB_position, UE_position] = channelModel(N, cellRadiusMax, cellRadiusMin, logNormalMean, logNormalDeviation)
% clc
% close all
% N: the number of mobile users 
% cellRadiusMax: maximum distance between a SeNB and a UE
% cellRadiusMin: mainimum distance between a SeNB and a UE

% % Sample parameters
% N = 6; 
% cellRadiusMax = 250; cellRadiusMin = 5;
% logNormalMean = 0; logNormalDeviation = 8;

% Coordinate of SeNBs
eNB_position = randsample(1:cellRadiusMax, 1); 
tera = 2*pi*rand(1, length(eNB_position));     % angle
xx_eNB = eNB_position.*cos(tera);             % x coordinate
yy_eNB = eNB_position.*sin(tera);             % y coordinate
BS_position_img = xx_eNB+1i.*yy_eNB;          % BS positions in the complex domain

% Coordinate of users - 5 random points
position_n1 = (1.5*(cellRadiusMax - cellRadiusMin).*rand(1,ceil(.5*N/10)) + 4*cellRadiusMin);
position_n2 = (1.5*(cellRadiusMax - cellRadiusMin).*rand(1,ceil(2.5*N/10)) + 4*cellRadiusMin);
position_n3 = (1.5*(cellRadiusMax - cellRadiusMin).*rand(1,ceil(3.5*N/10)) + 4*cellRadiusMin);
position_n4 = (1.5*(cellRadiusMax).*rand(1,ceil(1.5*N/10)) + 4*cellRadiusMin);
position_n5 = (1.5*(cellRadiusMax - cellRadiusMin).*rand(1,ceil(2*N/10)) + 4*cellRadiusMin);
% Coordinate of users - 5 random angles
theta1 = rand(1, length(position_n1));
theta2 = rand(1, length(position_n2)) + pi/2;
theta3 = (2*pi).*rand(1, length(position_n3)) + 5*pi/4;
theta4 = (2*pi).*rand(1, length(position_n4)) + 5*pi/4;
theta5 = (2*pi).*rand(1, length(position_n5)) + 0;
% Coordinate of users
[x1,y1] = pol2cart(theta1, position_n1);
[x2,y2] = pol2cart(theta2, position_n2);
[x3,y3] = pol2cart(theta3, position_n3);
[x4,y4] = pol2cart(theta4, position_n4);
[x5,y5] = pol2cart(theta5, position_n5);
UE_position = [x1+1i.*y1 x2+1i.*y2 x3+1i.*y3 x4+1i.*y4 x5+1i.*y5];
UE_position = UE_position(1:N);

% Sort users in the ascending order of the users' position
[~,IX] = sort(abs(UE_position));
UE_position = UE_position(IX);
xx_n = real(UE_position);
yy_n = imag(UE_position);
teta = angle(UE_position);

% Calculate the real distance between the eNodeB and mobile users
% dArray is a column vector
dArray = zeros(N,1);
for i = 1:N
    dArray(i) = sqrt((xx_n(i) - xx_eNB)^2 + (yy_n(i) - yy_eNB)^2);
end

% Caculate the channel gain between the eNodeB and mobile users 
% Channel gain hArray represents pathloss and log-normal shadowing
shadowing = zeros(N,1); 
for n = 1:N
    C = normrnd(logNormalMean, logNormalDeviation, 25, 1);
    shadowing(n) = mean(C);
end
% The path loss model for small cells is h = 140.7 + 36.7*log10(d)
hArray = -140.7 - 36.7.*log10(dArray/1000) + shadowing;
% Return the channel gain (normalized by the noise which is 95 dBm - 30 = 65 dB) in linear
hArray = db2lin(hArray - 10);

% %--------------Plotting the locations----------------%
%  LabledPlotting(xx_n, yy_n, teta,tera, UE_position, eNB_position, cellRadiusMin, cellRadiusMax)
end
