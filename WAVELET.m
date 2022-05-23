classdef WAVELET
    methods(Static)
        
        
        function [x,dt] = compute_signal(type_signal,dt)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%
        %%%%% Input:  - dt: the sampling time
        %%%%%         - type_signal: (0) gaussian pulse
        %%%%%                        (1) dirac
        %%%%%                        (2) noise
        %%%%%                        (3) cosine
        %%%%%                        (4) dirac+cosine
        %%%%%                        (5) parallel linear chirps
        %%%%%                        (6) sum of hyperbolic chirps
        %%%%%                        (7) signal Lewalle
        %%%%%                        (8) ramp
        %%%%%                        (9) diracs on the boundaries (coi)
        %%%%%
        %%%%% Output: - x : the signal 
        %%%%%         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        disp(' ')

        if(type_signal==0) % gaussian pulse
           disp('Signal: dirac')
           N=500;
           x=zeros(N,1);
           x=exp(-((1:N)-N/2).^2/50);

        elseif(type_signal==1) % dirac
           disp('Signal: dirac')
           N=200;
           dirac=zeros(N,1);
           dirac(N/2)=1/dt;
           x=dirac;

        elseif(type_signal==2) % noise
           disp('Signal: noise')
           x=randn(200,1);

        elseif(type_signal==3) % cosine
           disp('Signal: cosine')
           N=200;
           t=linspace(0,1,N);
           dt=diff(t);dt=dt(1);
           x=cos(2*pi*10*t);

        elseif(type_signal==4) % dirac + cosine
           disp('Signal: dirac+cosine')
           N=200;
           t=linspace(0,1,N);
           dt=diff(t);dt=dt(1);
           x=cos(2*pi*10*t);
           dirac=zeros(1,N);
           dirac(N/2)=1;
           x=x+dirac;

        elseif(type_signal==5) % parallel linear chirps
           disp('Signal: parralel linear chirps')
           N=1024;
           t=linspace(0,1,N);
           dt=diff(t);dt=dt(1);
           fs=1/dt;
           f1a=10;
           f2a=100;
           xa=exp(1i*2*pi*(f1a*t+(f2a-f1a)*t.^2/2));
           f1b=90;
           f2b=200;
           xb=exp(1i*2*pi*(f1b*t+(f2b-f1b)*t.^2/2));
           x=real(xa+xb);

        elseif(type_signal==6)
           disp('Signal: sum of hyperbolic chirps')
           N=2000;
           t=linspace(0,0.67,N);
           dt=diff(t);dt=dt(1);
           fs=1/dt;
           a1=1;
           a2=1;
           alpha1=150;
           alpha2=300;
           beta1=0.68;
           beta2=0.72;
           x=a1*cos(alpha1./(beta1-t))+a2*cos(alpha2./(beta2-t));

        elseif(type_signal==7) % signal Lewalle
           disp('Signal: Lewalle')
           load signal_lewalle
           x=y;dt=dx;
           clear dx y nx   

        elseif(type_signal==8) % ramp
           t=linspace(0,1,200);
           dt=diff(t);dt=dt(1);
           x=t;

        elseif(type_signal==9) % diracs on the boundaries --> cone of influence
           disp('Signal: dirac')
           N=200;
           dirac=zeros(N,1);
           dirac(1)=1;
           dirac(end)=1;
           x=dirac;   

        end

        disp(' ')
        disp('Subtract mean')
        x=x-mean(x);
        
        end
        
        
        function [SCALES] = compute_scale(smin,smax,nvoice)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%
        %%%%% Inputs: - smin: mininum scale.
        %%%%%         - smax: approximate maximum scale.
        %%%%%         - nvoice: number of sub-octaves per octave.
        %%%%%
        %%%%% Outputs: - SCALES: the vector containing scales for CWT
        %%%%%                    used as input for compute_cwt
        %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        I=fix((log2(smax/smin))*nvoice);
        SCALES = smin*2.^((0:I)/nvoice);
        
        end
        
        
        function [cwt_x] = compute_cwt(x,dt,SCALES,mother)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%
        %%%%% Continuous Wavelet Transform.
        %%%%%
        %%%%% Here we use the fact that: CWT(u,s)=(x*hs)(u) [where * is the convolution product]     
        %%%%% The CWT is calculated from: CWT= IFFT (FFT(x).FFT(hs)),
        %%%%%    where hs(t)= 1/sqrt(s).conj(psi(-t/s)) with psi the mother wavelet.
        %%%%%
        %%%%% The FFT of hs is obtained from the analytical FT after division by dt.   
        %%%%%
        %%%%% Inputs: - x: signal to be analysed, of size N.
        %%%%%         - dt: sampling time, of size N.
        %%%%%         - SCALES: the vector of scales at which the CWT is calculated, of size Nscale.     
        %%%%%         - mother: mother wavelet: either 'MORLET' or 'DOG'
        %%%%%
        %%%%%
        %%%%% Outputs - cwt_x: CWT of x. An array of size: Nscale x N,
        %%%%%                  where Nscale: number of scales (size of the array SCALES),
        %%%%%                  and   N: size of the signal.
        %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        N=length(x);
        time=0:dt:(N-1)*dt;
        x=reshape(x,1,N);
        Nscale=length(SCALES);



        %%%%% Fourier number
        f = WAVELET.compute_f(N,dt);


        %%%%% Wavelet Transform
        cwt_x = zeros(Nscale,N);
        xhat = fft(x);  
        for jscale = 1:Nscale   
            daughterhat=conj(WAVELET.compute_daughter_wavelet_hat(f,0,SCALES(jscale),mother));
            cwt_x(jscale,:) = ifft(xhat.*daughterhat);
        end
        cwt_x=cwt_x;
        
        end
        
        
        function [f] = compute_f(n,dt)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%   
        %%%%%   Returns an array with frequencies fk 
        %%%%%   The frequency vector corresponds to a n-point FFT
        %%%%%   with a sampling time dt.
        %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Compute the discrete frequencies
        f = [1:fix(n/2)];
        f = f.*(1/(n*dt));
        f = [0., f, -f(fix((n-1)/2):-1:1)];
        
        end
        

        function [daughter_hat] = compute_daughter_wavelet_hat(f,u,s,mother);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%
        %%%%%  Returns  - daughter_hat: the Fourier transform of the daughter wavelet
        %%%%%
        %%%%%  Inputs:  - f: frequency
        %%%%%           - u: time of analysis
        %%%%%           - s: scale
        %%%%%           - mother: 'Morlet'
        %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

        daughter_hat = sqrt(s)*WAVELET.compute_mother_wavelet_hat(s*f,mother).*exp(-1i*2*pi*f*u);
        
        end
        
        
        function [daughter] = compute_mother_wavelet_hat(f,mother);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%
        %%%%%
        %%%%%  Returns: - daughter: Fourier transform of the mother wavelet
        %%%%% 
        %%%%%
        %%%%%  Inputs:  - mother: 'Morlet' or 'DOG'
        %%%%%           - f: frequency
        %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        mother = upper(mother);
        omega=2*pi*f; % wavenumber

        if (strcmp(mother,'MORLET'))
            param=6.;
            omega0=param;
            norm=(4*pi)^0.25;
            daughter = norm*exp(-(omega-omega0).^2/2).*(omega>0);

        elseif (strcmp(mother,'DOG'))
            param=2;
            m=param;
            norm=(-1)^(m+1)/sqrt(gamma(m+0.5))*sqrt(2*pi);
            daughter=norm*((1i*omega).^m).*exp(-omega.^2/2);

        else
            error('Mother must be either MORLET or DOG!')
        end
        
        end
        
        

        function [xrecons]=compute_cwt_fftrecons(cwt_x,time,SCALES,Cpsi,mother)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%
        %%%%%    
        %%%%%  Return: - xrecons : reconstruction of signal from its Discrete-Time CWT.
        %%%%%
        %%%%%
        %%%%%  Inputs: - cwt_x : CWT of x, of size Nscale x N.
        %%%%%          - time  : time vector, of size N.
        %%%%%          - SCALES: vector of scales, of size Nscale.
        %%%%%          - Cpsi: admissibility condition.
        %%%%%          - mother: 'Morlet' or 'DOG'.
        %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        disp('Fast signal reconstruction')

        Nscale=length(SCALES);
        N=size(cwt_x,2);
        dt=time(2)-time(1);
        xrecons=zeros(1,N);
        dsvec=diff(SCALES); 

        % Frequency vector
        f = WAVELET.compute_f(N,dt);

        % Reconstruction
        for jscale=1:Nscale-1
           ds=dsvec(jscale);
           psis_hat=WAVELET.compute_mother_wavelet_hat(f*SCALES(jscale),mother);
           Ws_hat=fft(cwt_x(jscale,:));
           xrecons=xrecons+ifft(Ws_hat.*psis_hat)/(0.5*SCALES(jscale+1)+0.5*SCALES(jscale))^1.5*ds;
        end


        % Accounting for Cpsi
        if(strcmp(mother,'MORLET')) 
           xrecons=2/Cpsi*real(xrecons);
        elseif(strcmp(mother,'DOG')) 
           xrecons=1/Cpsi*real(xrecons); 
                                         
        end

        disp('Reconstruction completed')
        disp(' ')

        
        end
        


        function [xrecons]=compute_cwt_slowrecons(cwt_x,time,SCALES,Cpsi,mother)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%
        %%%%%
        %%%%%  Return: - xrecons: reconstruction of signal
        %%%%%
        %%%%%  Inputs: - cwt_x: CWT of size Nscale x N.
        %%%%%          - time: time vector, size N.
        %%%%%          - SCALES: scale vector, size Nscale.
        %%%%%          - Cpsi: constant characterizing the wavelet
        %%%%%          - mother: 'Morlet' or 'DOG'
        %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        disp('Signal reconstruction (long)...')

        Nscale=length(SCALES);
        N=size(cwt_x,2);
        dt=time(2)-time(1);

        xrecons=zeros(1,N);
        dsvec=diff(SCALES);

        for jscale=1:Nscale-1
           ds=dsvec(jscale);
           for jtime=1:N
              psi_us=WAVELET.compute_daughter_wavelet(time,time(jtime),SCALES(jscale),mother);
              xrecons=xrecons+psi_us*cwt_x(jscale,jtime)/(0.5*SCALES(jscale+1)+0.5*SCALES(jscale))^2*ds;
           end
        end
        xrecons=xrecons*dt;  

        if(strcmp(mother,'MORLET'))  
           xrecons=2/Cpsi*real(xrecons);  
        elseif(strcmp(mother,'DOG')) 
           xrecons=1/Cpsi*real(xrecons);       
        end

        disp('Reconstruction completed')
        disp(' ')
        
        end


        function f = compute_mother_wavelet(t,mother) 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%
        %%%%%   Returns: - f: mother wavelet at t.
        %%%%%
        %%%%%   Inputs: - t: running time.
        %%%%%           - mother: 'MORLET' or 'DOG'.
        %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        mother = upper(mother);         
        
        if (strcmp(mother,'MORLET'))
	        omega0=6;   % default value

            f=1/pi^(1/4)*exp(-t.^2/2).*exp(1i*omega0*t);
                
        elseif (strcmp(mother,'DOG'))
           m=2;   % default value
           
           COEFF=[1];ncoeff=1; % for m=0
           for i=1:m
              COEFF_dPndt=COEFF.*(0:1:ncoeff-1);
              COEFF_dPndt=[COEFF_dPndt(2:end) 0 0];
              
              ncoeff=ncoeff+1;
              COEFFnew=zeros(1,ncoeff);
              COEFFnew(2:end)=-COEFF; 
              COEFFnew=COEFFnew+COEFF_dPndt;
              COEFF=COEFFnew;
           end
           dogm=zeros(size(t));
           for i=1:length(COEFF)
              dogm=dogm+COEFF(i)*t.^(i-1); 
           end
           dogm=dogm.*exp(-t.^2/2);
             
               
           f=(-1)^(m+1)/sqrt(gamma(m+0.5))*dogm;
        
        else
	        error('Mother must be either MORLET or DOG!')
        end    
            
        end



        function daughter = compute_daughter_wavelet(t,u,s,mother) 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%
        %%%%%
        %%%%%   Returns: - daughter: daughter wavelet.
        %%%%%
        %%%%%   Inputs: - t: running time, in seconds
        %%%%%           - u: analysis time, in seconds
        %%%%%           - s: scale, in seconds
        %%%%%           - mother: 'MORLET' or 'DOG.'
        %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        daughter=1/sqrt(s)*WAVELET.compute_mother_wavelet((t-u)/s,mother);
        
        end
        
    end
    
end