% Script based on original MAIN.m

addpath("..")
clear all
%% read info data (two examples are given here)
% info_test4;
info_test5;

load(B1_file); % load B1 map
Nch = size(B1,2); 
frequencies = [ones(1,singles), const*ones(1,intervals)];
figure;bar(frequencies);title('tip angles repetitions (frequencies)');grid on;
Nt = sum(frequencies);% number of time points
klim = min([klim, round(Nt/2)]);

%% set angles for target: 
alpha0 = transpose(a0(:));
c = zeros(1,Nt+1);
c(initial:end) = 1; 
[F] = EPG_forward1(alpha0,'ESP',ESP,'T1',T1,'T2',T2); % set target, first simulate for ideal alpha
target = full(F(2,:)); % target F vector (containes Nt-1 echoes)
target = abs(target);
Ns = size(B1,1);


%% starting value. Start with quadrature values
alpha_start = kron(ones(1,Nch),alpha0(cumsum(frequencies)));
phi_start = kron(zeros(1,Nch),ones(1,length(alpha0(cumsum(frequencies)))));
TARGET = repmat(target(:),[ 1 Ns]); % target for all spatial locations
param_start = [real(alpha_start.*exp(1i*phi_start)), imag(alpha_start.*exp(1i*phi_start))];

%%
param_start1 = red2full(param_start,frequencies,Nch);
[obj, grad,FF] = obj_EPG11(param_start1,ESP,T1,T2,c,B1,TARGET);
N = size(FF,1)/2;
param_start1 = red2full(param_start,frequencies,Nch);
xstart = reshape(param_start1(1:Nt*Nch),Nt,Nch).'; %real parts
ystart = reshape(param_start1(Nt*Nch+1:end),Nt,Nch).'; % imaginary parts
Alpha_start = abs(xstart+1i*ystart);
Phi_start = angle(xstart+1i*ystart);

%% use the new thetaeff

thetain = zeros(Ns,Nt,Nch);
for ch = 1:Nch
    thetain(:,:,ch) = B1(:,ch)*(xstart(ch,:)+1i*ystart(ch,:)); % effective theta
end
thetain(:,1,:) = pi/16*ones(Ns,1,Nch);
[obj1, grad1,FF1] = obj_EPG11_fromeff(thetain,ESP,T1,T2,c,B1,TARGET);

%%
figure(91)
subplot(5,2,1);
plot(180/pi*Alpha_start.','LineWidth',2);ylim([0 180*max_a]);
title('STARTING \alpha');
subplot(5,2,3);
plot(180/pi*Phi_start.','LineWidth',2);ylim([-180 180]);
title('STARTING \phi');
subplot(5,2,5);
plot(squeeze(abs((FF(2,:,:)+1i*FF(N+2,:,:)))),'LineWidth',.5);title('STARTING SIGNAL');xlim([1 Nt+1]);grid on;
hold on;plot(abs(sum(TARGET,2)/Ns),'bo','LineWidth',2);
subplot(5,2,7);
imagesc(abs(squeeze(FF(2,:,:)+1i*FF(N+2,:,:))).');title('STARTING SIGNAL');grid on;ylabel('voxels')
subplot(5,2,9);
imagesc(angle(squeeze(FF(2,:,:)+1i*FF(N+2,:,:))).');title('STARTING SIGNAL');grid on;ylabel('voxels')

subplot(5,2,2);
plot(180/pi*Alpha_start.','LineWidth',2);ylim([0 180*max_a]);
title('STARTING \alpha');
subplot(5,2,4);
plot(180/pi*Phi_start.','LineWidth',2);ylim([-180 180]);
title('STARTING \phi');
subplot(5,2,6);
plot(squeeze(abs((FF1(2,:,:)+1i*FF1(N+2,:,:)))),'LineWidth',.5);title('STARTING SIGNAL');xlim([1 Nt+1]);grid on;
hold on;plot(abs(sum(TARGET,2)/Ns),'bo','LineWidth',2);
subplot(5,2,8);
imagesc(abs(squeeze(FF1(2,:,:)+1i*FF1(N+2,:,:))).');title('STARTING SIGNAL');grid on;ylabel('voxels')
subplot(5,2,10);
imagesc(angle(squeeze(FF1(2,:,:)+1i*FF1(N+2,:,:))).');title('STARTING SIGNAL');grid on;ylabel('voxels')


