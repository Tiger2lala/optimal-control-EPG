% Test the importing .mat pulse file. 
% Script based on original MAIN.m

addpath("..")
addpath("../../BlochSimDYBZMH/")
%% info
clear all
ESP = 12.0 ;%<- echo spacing, ms
T1 = 1500; %<- T1,ms
T2 = 100;  %<- T2,ms
klim = 25; % order of maximum k, use 25
max_a = 1.5; % maximum flip = max_a*pi;
red_b1 = 1; % spatial undersampling factor for B1 maps

singles = 2; % number of individual pulses at the beginning of sequence
const = 10; % interval size of constant pulses after the individual ones
intervals = 0; % number of intervals with constant flips
initial = 3;

Adjfile = "../../volunteerAdj/AdjDataUser.mat";
load(Adjfile, 'Adj')
ROI = zeros(size(Adj.W));
ROI(:,:,36) = Adj.W(:,:,36);
ROI = logical(ROI);
B1 = b1convert(Adj, 245, ROI);

a0 = pi/180*[90 90]; 

%%
Nch = size(B1,2); 
frequencies = [ones(1,singles), const*ones(1,intervals)];
Nt = sum(frequencies);% number of time points
klim = min([klim, round(Nt/2)]);
alpha0 = transpose(a0(:));
c = zeros(1,Nt+1);
c(initial:end) = 1; 
[F] = EPG_forward1(alpha0,'ESP',ESP,'T1',T1,'T2',T2); % set target, first simulate for ideal alpha
target = full(F(2,:)); % target F vector (containes Nt-1 echoes)
target = abs(target);
Ns = size(B1,1);

alpha_start = kron(ones(1,Nch),alpha0(cumsum(frequencies)));
phi_start = kron(zeros(1,Nch),ones(1,length(alpha0(cumsum(frequencies)))));
TARGET = repmat(target(:),[ 1 Ns]); % target for all spatial locations
param_start = [real(alpha_start.*exp(1i*phi_start)), imag(alpha_start.*exp(1i*phi_start))];

%%
param_start1 = red2full(param_start,frequencies,Nch);
xstart = reshape(param_start1(1:Nt*Nch),Nt,Nch).'; %real parts
ystart = reshape(param_start1(Nt*Nch+1:end),Nt,Nch).'; % imaginary parts

thetain = zeros(Ns,Nt,Nch);
for ch = 1:Nch
    thetain(:,:,ch) = B1(:,ch)*(xstart(ch,:)+1i*ystart(ch,:)); % effective theta
end
pulsefile = "../../ManualPulse/CP-1.5mm-exc.mat";
thetain(:,1,:) = load_pulse_theta(Adjfile, pulsefile, ROI);
[obj1, grad1,FF1] = obj_EPG11_fromeff(thetain,ESP,T1,T2,c,B1,TARGET);

%%
figure
subplot(2,2,1)
outimg = zeros(55,55);
outimg(ROI(:,:,36))=asin(abs(FF1(2,2,:)+1i*FF1(5,2,:)));
imagesc(outimg); axis equal tight
title('CP mode acrsin(EPG signal) after 90 exc')

subplot(2,2,3)
outimg(ROI(:,:,36))=(angle(FF1(2,2,:)+1i*FF1(5,2,:)));
imagesc(outimg); axis equal tight; clim([-pi pi])

subplot(2,2,2)
simimg = zeros(55,55);
simimg(ROI(:,:,36))=abs(thetain(:,1,1));
imagesc(simimg); axis equal tight
title('Bloch sim FA')

subplot(2,2,4)
simimg(ROI(:,:,36))=angle(thetain(:,1,1));
imagesc(simimg); axis equal tight
