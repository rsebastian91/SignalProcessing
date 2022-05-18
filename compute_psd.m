%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%
%%%%%
%%%%% Calculate Power Spectral Density using 
%%%%% the raw and averaged periodogram
%%%%% for cosine and white noise
%%%%%
%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%%%%%% I. SIGNAL DEFINITION

%%% GENERAL
nbdata=8192;                     %%% number of data points
Ts=0.001;                        %%% sampling period
fs=1/Ts;                         %%% sampling frequency
tmax=(nbdata-1)*Ts;              %%% total duration of the signal

%%% COSINE PART
freq=100;                        %%% frequency of the cosine
period=1/freq;                   %%% period of the cosine
t=0:Ts:tmax ;                    %%% time vector
A=0;                             %%% amplitude of the cosine     
x=A*cos(2*pi*freq*t);            %%% cosine

%%% ADDITION OF NOISE
mean_noise=0;
std_noise=5;                     %%% standard deviation (sigma) 

noise=mean_noise+std_noise*randn(size(x));
x=x+noise;

%%%%% FIGURE 1 - SIGNAL
fig = figure(1);
orient portrait;

set(fig,'papertype','a4letter');
set(fig,'units','centimeters','paperunits','centimeters');
set(fig,'paperposition',[2 8 14 10]);
set(fig,'position',[2 8 14 10]);

hold on
plot(t(1:100),x(1:100),'b')
xlabel('t (s)','fontsize',16)
ylabel('x(t)','fontsize',16)



%%%%% II. PSD CALCULATED WITH THE RAW PERIODOGRAM
% Number of points used for the raw periodogram (1<Nraw<nbdata)
Nraw=nbdata;

% Windowing function
iwindow=1;  % 1 - Rectangular, 2 - Hanning, 3 - Blackman

[f_raw,Sxx_raw]=PSD.compute_raw_psd(x(1:Nraw),1/Ts,iwindow);

disp(' ')
disp(['Mean of raw periodogram: ' num2str(mean(Sxx_raw))])
disp(['Standard deviation of raw periodogram: ' num2str(std(Sxx_raw))])
disp(' ')

%%%%%% III. PSD CALCULATED WITH THE AVERAGED PERIODOGRAM
% Number of points per block for the Welch periodogram (1<Nwelch<nbdata)
Nwelch=300;                    
[f_welch,Sxx_welch,Mwelch]=PSD.compute_welch(x,1/Ts,Nwelch,iwindow);

disp(' ')
disp(['Mean of averaged periodogram: ' num2str(mean(Sxx_welch))])
disp(['Standard deviation of averaged periodogram: ' num2str(std(Sxx_welch))])
disp(' ')

%%%%% FIGURE 2 - PSD
fig = figure(2);
orient portrait;

set(fig,'papertype','a4letter');
set(fig,'units','centimeters','paperunits','centimeters');
set(fig,'paperposition',[3 7 14 10]);
set(fig,'position',[3 7 14 10]);

hold on
xlabel('f (Hz)','fontsize',16)
ylabel('S_{xx}(f)','fontsize',16)
plot(f_raw,Sxx_raw,'b');
plot(f_welch,Sxx_welch,'g');
legend('Raw periodogram','Averaged periodogram')


%%%%% IV. ENERGY CALCULATION
% Time domain 
Pt=PSD.compute_power_time(x);
disp(['Power calculated in the time domain: Pt=' num2str(Pt)]);
% Frequency domain / raw
Pf_raw=PSD.compute_power_freq(Sxx_raw);
disp(['Power calculated in the frequency domain: Pf_raw=' num2str(Pf_raw)]);
% Frequency domain / averaged
Pf_welch=PSD.compute_power_freq(Sxx_welch);
disp(['Power calculated in the frequency domain: Pf_welch=' num2str(Pf_welch)]);
disp(' ')

%%%%% CHECK CONVERGENCE
% The vector Nvec contains several block sizes
Nvec=[100:200:4000];   
% For each of these block size, we find the corresponding number
% of blocks, the mean of Sxx_welch, and the standard deviation of
% Sxx_welch.
[Mvec,mean_Sxx,sigma_Sxx]= PSD.compute_error(x,fs,Nvec,iwindow);

%%%%% FIGURE 3 - PSD ESTIMATE MEAN
fig = figure(3);
orient portrait;

set(fig,'papertype','a4letter');
set(fig,'units','centimeters','paperunits','centimeters');
set(fig,'paperposition',[4 6 14 10]);
set(fig,'position',[4 6 14 10]);

hold on
xlabel('M (number of blocks)','fontsize',16)
ylabel('Mean [ S_{xx} ]','fontsize',16)
plot(Mvec,mean_Sxx,'b');
yline(std_noise^2,'g')



%%%%% FIGURE 4 - PSD ESTIMATE STANDARD VARIATION
fig = figure(4);
orient portrait;

set(fig,'papertype','a4letter');
set(fig,'units','centimeters','paperunits','centimeters');
set(fig,'paperposition',[5 5 14 10]);
set(fig,'position',[5 5 14 10]);

hold on
xlabel('M (number of blocks)','fontsize',16)
ylabel('\sigma [ S_{xx} ]','fontsize',16)
plot(Mvec,sigma_Sxx,'b');
plot((1:200),std_noise^2./sqrt((1:200)),'g')

leg=legend('$\sigma[S_{xx}]$ (observed)','$\sigma^2/\sqrt{M}$ (theory)');
set(leg,'interpreter','latex')



