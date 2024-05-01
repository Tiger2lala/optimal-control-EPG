% Test the use obj_EPG13 wrapper in optimizer. 
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
pow_factor = 2.0; % scale total power factor
chpow_factor = 1.5;
%%
singles = 5; % number of individual pulses at the beginning of sequence
const = 3; % interval size of constant pulses after the individual ones
intervals = 0; % number of intervals with constant flips
initial = 3;

Adjfile = "../../phantomAdj/AdjDataUser.mat";
load(Adjfile, 'Adj')
ROI = zeros(size(Adj.W));
ROI= Adj.W;
ROI = logical(ROI);
B1 = b1convert(Adj, 245, ROI);

a0 = pi/180*[90 180 180 180 180]; 

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

%% prepare optimize
maxiter = 20;
param_all = red2full(param_start,frequencies,Nch);
max_power = pow_factor*norm(param_all(:))^2; % constraint power to power of starting sequence
max_chpow = chpow_factor*ones(1,Nch)*sum(alpha0.^2);
[obj, grad,FF] = obj_EPG13(param_start,ESP,T1,T2,c,B1,TARGET,frequencies,10e+06);
N = size(FF,1)/2;
options = optimoptions('fmincon','MaxIter',maxiter, 'MaxFunEvals',1000*length(alpha_start)*5, ...
    'GradObj','on','GradConstr','on','DerivativeCheck','off','TolFun',1.0e-09); 
options = optimoptions(options,'Display','iter-detailed');
options = optimoptions(options,'PlotFcns',@optimplotfval);
Z_TARGET = squeeze(angle(FF(2,:,:)+1i*FF(N+2,:,:)));
Z_TARGETs = unwrap(Z_TARGET,[],1);
Z_TARGETs(2:end,:) = repmat(mean(Z_TARGETs(initial:end,:),1),[Nt 1]);
TARGET = abs(TARGET).*exp(1i*Z_TARGETs);

%% optimize original
[sol1 fval] = fmincon(@(thet) obj_EPG13(thet,ESP,T1,T2,c,B1,TARGET,frequencies,klim),param_start,[],[],[],[],[],[],@(thet) limit_RF(thet,max_a*pi,Nch),options); 

sol = red2full(sol1,frequencies,Nch);
param_start1 = red2full(param_start,frequencies,Nch);
xopt = reshape(sol(1:Nt*Nch),Nt,Nch).'; %real parts
yopt = reshape(sol(Nt*Nch+1:end),Nt,Nch).'; % imaginary parts
[obj, grad,FF1] = obj_EPG11(sol,ESP,T1,T2,c,B1,TARGET);
state1 = squeeze(FF1(2,:,:)+1i*FF1(N+2,:,:));
alpha_opt = abs(xopt+1i*yopt);
phi_opt = angle(xopt+1i*yopt);

%% from read-in excitation
param_start_red = reshape(param_start, [],Nch,2);
param_start_red = param_start_red(2:end,:,:);
pulsefile = "../../ManualPulse/CP-1.5mm-exc.mat";
exc = load_pulse_theta(Adjfile, pulsefile, ROI);

[sol2 fval] = fmincon(@(thet) obj_EPG13_exc_wrapper(thet,exc,ESP,T1,T2,c,B1,TARGET,frequencies,klim),param_start_red(:),[],[],[],[],[],[],@(thet) limit_RF(thet,max_a*pi,Nch),options); 
sol2 = reshape(sol2, [],Nch,2);
sol2(2:end+1, :, :) = sol2;
sol2 = sol2(:);
sol = red2full(sol2,frequencies,Nch);
param_start1 = red2full(param_start,frequencies,Nch);
xopt = reshape(sol(1:Nt*Nch),Nt,Nch).'; %real parts
yopt = reshape(sol(Nt*Nch+1:end),Nt,Nch).'; % imaginary parts

thetain = zeros(Ns,Nt,Nch);
for ch = 1:Nch
    thetain(:,:,ch) = B1(:,ch)*(xopt(ch,:)+1i*yopt(ch,:)); % effective theta
end
thetain(:,1,:) = exc;
[obj, grad,FF2] = obj_EPG11_fromeff(thetain,ESP,T1,T2,c,B1,TARGET);
state2 = squeeze(FF1(2,:,:)+1i*FF1(N+2,:,:));
alpha_opt2 = abs(xopt+1i*yopt);
phi_opt2 = angle(xopt+1i*yopt);

%% starting
xstart = reshape(param_start1(1:Nt*Nch),Nt,Nch).'; %real parts
ystart = reshape(param_start1(Nt*Nch+1:end),Nt,Nch).'; % imaginary parts
[obj, grad,FF] = obj_EPG11(param_start1,ESP,T1,T2,c,B1,TARGET);
Alpha_start = abs(xstart+1i*ystart);
Phi_start = angle(xstart+1i*ystart);

%%


figure
subplot(5,3,1);
plot(180/pi*Alpha_start.','LineWidth',2);ylim([0 180*max_a]);
title('STARTING \alpha');
subplot(5,3,4);
plot(180/pi*Phi_start.','LineWidth',2);ylim([-180 180]);
title('STARTING \phi');
subplot(5,3,7);
plot(squeeze(abs((FF(2,:,:)+1i*FF(N+2,:,:)))),'LineWidth',.5);title('STARTING SIGNAL');xlim([1 Nt+1]);grid on;
hold on;plot(abs(sum(TARGET,2)/Ns),'bo','LineWidth',2);

subplot(5,3,2);
plot(180/pi*alpha_opt.','LineWidth',2);ylim([0 180*max_a]);
title('Default OPTIMIZED \alpha');
subplot(5,3,5);
plot(180/pi*phi_opt.','LineWidth',2);ylim([-180 180]);
title('Default OPTIMIZED \phi');
subplot(5,3,8);
plot(squeeze(abs((FF1(2,:,:)+1i*FF1(N+2,:,:)))),'LineWidth',.5);title('default OPTIMIZED SIGNAL');xlim([1 Nt+1]);grid on;
hold on;plot(abs(sum(TARGET,2)/Ns),'bo','LineWidth',2);

subplot(5,3,3);
plot(180/pi*alpha_opt2.','LineWidth',2);ylim([0 180*max_a]);
title('Exc OPTIMIZED \alpha');
subplot(5,3,6);
plot(180/pi*phi_opt2.','LineWidth',2);ylim([-180 180]);
title('Exc OPTIMIZED \phi');
subplot(5,3,9);
plot(squeeze(abs((FF2(2,:,:)+1i*FF2(N+2,:,:)))),'LineWidth',.5);title('Exc OPTIMIZED SIGNAL');xlim([1 Nt+1]);grid on;
hold on;plot(abs(sum(TARGET,2)/Ns),'bo','LineWidth',2);


subplot(5,3,10);
imagesc(abs(squeeze(FF(2,:,:)+1i*FF(N+2,:,:))).');title('STARTING SIGNAL');grid on;ylabel('voxels')
subplot(5,3,11);
imagesc(abs(squeeze(FF1(2,:,:)+1i*FF1(N+2,:,:))).');title('OPTIMIZED SIGNAL');grid on;ylabel('voxels')
subplot(5,3,12);
imagesc(abs(squeeze(FF2(2,:,:)+1i*FF2(N+2,:,:))).');title('EXC OPTIMIZED SIGNAL');grid on;ylabel('voxels')
subplot(5,3,13);
imagesc(angle(squeeze(FF(2,:,:)+1i*FF(N+2,:,:))).');title('STARTING SIGNAL');grid on;ylabel('voxels')
subplot(5,3,14);
imagesc(angle(squeeze(FF1(2,:,:)+1i*FF1(N+2,:,:))).');title('OPTIMIZED SIGNAL');grid on;ylabel('voxels')
subplot(5,3,15);
imagesc(angle(squeeze(FF2(2,:,:)+1i*FF2(N+2,:,:))).');title('EXC OPTIMIZED SIGNAL');grid on;ylabel('voxels')

%%
img = zeros(Adj.image_m, Adj.image_n);

% signal at half the train?
figure
subplot(1,3,1)
img(ROI) = abs(squeeze(FF(2,ceil(end/2),:)+1i*FF(N+2,ceil(end/2),:))).';
imagesc(img);title("Start halfway signal");axis equal tight

subplot(1,3,2)
img(ROI) = abs(squeeze(FF1(2,ceil(end/2),:)+1i*FF1(N+2,ceil(end/2),:))).';
imagesc(img);title("Optimized halfway signal");axis equal tight

subplot(1,3,3)
img(ROI) = abs(squeeze(FF2(2,ceil(end/2),:)+1i*FF2(N+2,ceil(end/2),:))).';
imagesc(img);title("Exc+optimized halfway signal");axis equal tight
