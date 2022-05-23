%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%
%%%%%
%%%%%  Discrete-Time Continuous Wavelet Transform (CWT)
%%%%%
%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

nvoice=20;                         % number of sub-oactves per octave
mother='DOG';                      % 'Morlet' or 'DOG' 
o_slow_reconstruction_wavelet=0;   

% Load signal
type_signal=1;  % between 0 and 9
dt=0.001;
x=WAVELET.compute_signal(type_signal,dt);
N=length(x);
time=0:dt:(N-1)*dt;

% Display parameter setting
disp(' ')
disp(['N=' num2str(N) '    dt=' num2str(dt)])

mother=upper(mother);
disp(' ')
disp(['Mother: ' mother])

disp(' ')
if(o_slow_reconstruction_wavelet==1)
   disp('Slow reconstruction using wavelets')
else
   disp('Fast reconstruction using wavelets')
end



% Constant Cpsi (admissibility condition)
mother=upper(mother);
if (strcmp(mother,'MORLET'))
   Cpsi=1.06;
   Te=0.707;
elseif (strcmp(mother,'DOG'))
   Cpsi=2.36;
   Te=1.08;
else
   disp('STOP: Cpsi not defined!')
end

disp(' ')
disp(['Cpsi=' num2str(Cpsi)]);
disp(' ')


%%%%% SCALES
smin=0.1*dt;    
smax=3*N*dt;    
SCALES=WAVELET.compute_scale(smin,smax,nvoice);


%%%%% CWT
[cwt_x]=WAVELET.compute_cwt(x,dt,SCALES,mother); 

%%%%% SCALOGRAM
power = (abs(cwt_x)).^2 ;


disp('Full reconstruction')
%%%%% FAST RECONSTRUCTION
xr=WAVELET.compute_cwt_fftrecons(cwt_x,time,SCALES,Cpsi,mother);

%%%%% SLOW RECONSTRUCTION
if(o_slow_reconstruction_wavelet)  
   xr2=WAVELET.compute_cwt_slowrecons(cwt_x,time,SCALES,Cpsi,mother);
end


disp('Large scale reconstruction')
%%%%% SMALL SCALE REMOVAL
cwt_x(1:120,:)=0;
xls=WAVELET.compute_cwt_fftrecons(cwt_x,time,SCALES,Cpsi,mother);


%%%%% FIGURES
startx=2;
starty=6;
sizex=14;
sizey=11;

%%%%% FIGURE 1 - REAL PART OF CWT - LINEAR SCALE
fig = figure(1);
orient portrait;
set(fig,'papertype','a4letter');
set(fig,'units','centimeters','paperunits','centimeters');
set(fig,'paperposition',[startx starty sizex sizey]);
set(fig,'position',[startx starty sizex sizey]);startx=startx+0.5;starty=starty-0.5;

hold on
pcolor(time,SCALES,real(cwt_x))
set(gca, 'FontSize', 12);
colormap(jet)
shading flat
xlim([time(1) time(end)])
ylim([SCALES(1) SCALES(end)])
xlabel('t (s)','fontsize',16)
ylabel('scale (s)','fontsize',16)
title(['REAL PART of the CWT of x(t)'],'color','b','fontweight','bold')
box on




%%%%% FIGURE 2- REAL PART OF CWT - LOG SCALE
fig = figure(2);
orient portrait;
set(fig,'papertype','a4letter');
set(fig,'units','centimeters','paperunits','centimeters');
set(fig,'paperposition',[startx starty sizex sizey]);
set(fig,'position',[startx starty sizex sizey]);startx=startx+0.5;starty=starty-0.5;

pcolor(time,log2(SCALES/dt),real(cwt_x))
set(gca, 'FontSize', 12);
colormap(jet)
shading flat
hold on
%%% Cone of influence
plot(time,log2(min(time/3/Te,(time(end)-time)/3/Te)/dt),'k');
xlim([time(1) time(end)])
xlabel('t (s)','fontsize',16)
ylabel('log_2(scale/dt)','fontsize',16)
title(['REAL PART of the CWT of x(t)'],'color','b','fontweight','bold')
box on



%%%%% FIGURE 3 - SCALOGRAM=ABS(CWT)^2 - LINEAR SCALE
fig = figure(3);
orient portrait;
set(fig,'papertype','a4letter');
set(fig,'units','centimeters','paperunits','centimeters');
set(fig,'paperposition',[startx starty sizex sizey]);
set(fig,'position',[startx starty sizex sizey]);startx=startx+0.5;starty=starty-0.5;

pcolor(time,log2(SCALES/dt),power)
set(gca, 'FontSize', 12);
colormap(jet)
shading flat
xlim([time(1) time(end)])
xlabel('t (s)','fontsize',16)
ylabel('log_2(scale/dt)','fontsize',16)
title(['SCALOGRAM |CWT|^2'],'color','b','fontweight','bold')
box on


%%%%% FIGURE 4 - RECONSTRUCTION
fig = figure(4);
orient portrait;
set(fig,'papertype','a4letter');
set(fig,'units','centimeters','paperunits','centimeters');
set(fig,'paperposition',[startx starty sizex sizey]);
set(fig,'position',[startx starty sizex sizey]);startx=startx+0.5;starty=starty-0.5;

%
hold on
plot(time,x,'r')
plot(time,xr,'+b')
plot(time,xls,'k','linewidth',2)
if(o_slow_reconstruction_wavelet) 
   plot(time,xr2,':r')
end
set(gca, 'FontSize', 10);
xlabel('t (s)','fontsize',16)
ylabel('x, x_r','fontsize',16)
legend('original','fast wavelet reconstruction','large scale recons','slow wavelet reconstruction')
title('Reconstruction','color','b','fontweight','bold')



%%%%% FIGURE 5 - CHECKING WHEN SIGNAL = DIRAC
if(type_signal==1)
   fig = figure(5);
   orient portrait;
   set(fig,'papertype','a4letter');
   set(fig,'units','centimeters','paperunits','centimeters');
   set(fig,'paperposition',[startx starty sizex sizey]);
   set(fig,'position',[startx starty sizex sizey]);startx=startx+0.5;starty=starty-0.5;
   
   % real part
   subplot(2,1,1)
   hold on
   I=find(x>0.5,1);
   tdirac=time(I);                   %%% time of the dirac
   scale_plot=10*dt;                 %%% scale at which the comparison is done
   iscale= find(SCALES>scale_plot,1); %%% index of the scale array corresponding to that scale
   cwt_anly=conj(WAVELET.compute_daughter_wavelet(-time,-tdirac,SCALES(iscale),mother));
   xlabel('t (s)','fontsize',16)
   ylabel('Real(CWT)','fontsize',16)
   title(['Comparison CWT / analysing wavelet   -    scale=' num2str(SCALES(iscale))],'color','b')
   hold on
   plot(time,real(cwt_x(iscale,:)),'*k') % dummy plot for legend not to return error when no plot
   plot(time,real(cwt_anly),'-r','LineWidth',2)
   set(gca, 'FontSize', 10);
   legend('CWT','Analysing wavelet')
   % imaginary part
   subplot(2,1,2)
   
   hold on
   xlabel('t (s)','fontsize',16)
   ylabel('Imag(CWT)','fontsize',16)
   title(['Comparison CWT / analysing wavelet   -    scale=' num2str(SCALES(iscale))],'color','b')
   plot(time,imag(cwt_x(iscale,:)),'*k') % dummy plot for legend not to return error when no plot
   plot(time,imag(cwt_anly),'-r','LineWidth',2)
   set(gca, 'FontSize', 10);
   legend('CWT','Analysing wavelet')
end
