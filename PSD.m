classdef PSD
    methods(Static)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%
        %%%%%
        %%%%% Returns the Power Spectral Density using the raw periodogram
        %%%%%
        %%%%%
        %%%%% Inputs: - x  : input signal 
        %%%%%         - fs : sampling frequency
        %%%%%         - iwindow: type of window used: (1):  Rectangular
        %%%%%                                         (2):  Hanning
        %%%%%                                         (3):  Blackman
        %%%%%
        %%%%%
        %%%%% Outputs:  - f   : a frequency vector spanning [-fs/2 fs/2[
        %%%%%           - Sxx : the PSD vector
        %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        function [f,Sxx] = compute_raw_psd(x,fs,iwindow)


        %%%%% SIZE OF THE SIGNAL
        N=length(x);

        %%%%% WINDOWING FUNCTION
        disp(' ')
        if(iwindow==1)
           w=ones(N,1)';
           disp('WINDOW: rectangular')
        elseif(iwindow==2)
           w=hanning(N)';
           disp('WINDOW: Hanning')
        elseif(iwindow==3)
           w=blackman(N)';   
           disp('Window: Blackman')
        end
        Cw=N/sum(w.^2);   %%% Corrective factor

        %%%%% FREQUENCY VECTOR
        f=((-N/2):((N/2)-1))*fs/N;

        %%%%% POWER SPECTRAL DENSITY
        Sxx=(Cw*(abs(fftshift(fft(x.*w)))).^2)/N;
        
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%
        %%%%%
        %%%%% Power Spectral Density using Welch's method.
        %%%%% A 50% overlap is used between blocks.
        %%%%%
        %%%%% Inputs: - x  : input signal 
        %%%%%         - fs : sampling frequency
        %%%%%         - N  : size of a block
        %%%%%         - iwindow: type of window used: (1):  Rectangular
        %%%%%                                         (2):  Hanning
        %%%%%                                         (3):  Blackman
        %%%%%
        %%%%%
        %%%%% Outputs:  - f   : a frequency vector going from -fs/2 to fs/2
        %%%%%           - Sxx : the PSD vector
        %%%%%           - M   : the number of blocks 
        %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [f,Sxx,M] = compute_welch(x,fs,N,iwindow)


        %%%%% CHECKING BLOCK SIZE
        if(N>length(x))
           disp(['*** ERROR: Block size N should not be larger than the signal size.'])
           disp('*** STOP !!!')
           return
        end


        %%%%% WINDOWING FUNCTION
        if(iwindow==1)
           w=ones(N,1)';
           %disp('WINDOW: rectangular')
        elseif(iwindow==2)
           w=hanning(N)';
           %disp('WINDOW: Hanning')
        elseif(iwindow==3)
           w=blackman(N)';   
           %disp('Window: Blackman')
        end
        Cw=N/sum(w.^2);   %%% Corrective factor


        %%%%% CALCULATION OF THE NUMBER OF BLOCKS (M)
        noverlap=floor(N/2);   % 50% overlap
        go_on=1;
        n_end_block=N;
        M=1;
        while(go_on==1)
           n_end_block=n_end_block+noverlap;
           if(n_end_block<=length(x))
              M=M+1;
              go_on=1;
           elseif(n_end_block>length(x))
              go_on=0;
           end
        end
        disp(['Welch: number of blocks: M= ' num2str(M)]);


        %%%%% FREQUENCY VECTOR
        f=((-N/2):((N/2)-1))*fs/N;

        %%%%% CALCULATION OF THE AVERAGED PSD (Sxx)
        %%%%% WITH THE PSD OF EACH BLOCK (Sxx_current_block)      
        Sxx=zeros(1,N);

        for j=1:M
         Sxx=Sxx+(Cw*abs(fftshift(fft(x(1+(j-1)*noverlap:N+(j-1)*noverlap).*w))).^2)/N;   
        end
        Sxx=Sxx/M;
        
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%
        %%%%%
        %%%%% Returns power of the temporal signal
        %%%%%
        %%%%% Inputs: -    x  : input signal 
        %%%%%
        %%%%% Outputs:  - P   : signal power
        %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function P = compute_power_time(x)

        N=length(x);
        P=sum(x.^2)/N;
        
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%
        %%%%%
        %%%%% Returns power of the PSD
        %%%%%
        %%%%% Inputs: -  Sxx  : PSD 
        %%%%%
        %%%%% Outputs:  - P   : Power
        %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function P = compute_power_freq(Sxx)


        N=length(Sxx);
        P=sum(Sxx)/N;
        
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%
        %%%%%
        %%%%% Calulate the standard deviation of the PSD estimate.
        %%%%%
        %%%%% Uses: - compute_welch
        %%%%%
        %%%%% Inputs: - x       : input signal 
        %%%%%         - fs      : sampling frequency
        %%%%%         - Nvec    : a vector containing different block sizes
        %%%%%         - iwindow : type of window used: (1):  Rectangular
        %%%%%                                          (2):  Hanning
        %%%%%                                          (3):  Blackman
        %%%%%
        %%%%% Outputs:  - Mvec       : a vector of size Nvec containing the number of blocks used
        %%%%%           - mean_Sxx   : a vector of size Nvec containing the mean of the PSD
        %%%%%           - std_Sxx    : a vector of size Nvec containing the std deviation of the PSD  
        %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



        function [Mvec,mean_Sxx,std_Sxx]= compute_error(x,fs,Nvec,iwindow)


        Mvec=zeros(1,length(Nvec));
        mean_Sxx=zeros(1,length(Nvec));
        std_Sxx=zeros(1,length(Nvec));


        for i=1:length(Nvec)

            [f_welch,Sxx_welch,Mvec(i)]=PSD.compute_welch(x,fs,Nvec(i),iwindow);
            mean_Sxx(i)=mean(Sxx_welch);
            std_Sxx(i)=std(Sxx_welch);
        end
        
        end
        
    end
    
end