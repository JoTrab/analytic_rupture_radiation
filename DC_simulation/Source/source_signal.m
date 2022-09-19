function [signal, deriv] = source_signal(Source,Params)
%% same function for both subgits
%%
type_switch = Source.type_switch ;
figure_switch = Source.figure_switch ;
time = Params.time ;
delta_t = Params.delta_t ;
Mo = Source.Magnitude;    %scalar moment
switch(type_switch)
    case{'sin_mod'}
        clear Ntime signal frequence spectrum
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Définition du signal d'émission
        Fc=Source.Frequency_Hz ;
        Fmax=Fc + ceil(0.5*Source.Frequency_Hz*Source.Bandwidth);
        Fmin=Fc - ceil(0.5*Source.Frequency_Hz*Source.Bandwidth);
        omc=2*pi*Fc;
        
        % Tmax=0.1;
        
        Ntime=length(time);
        % if rem(Ntime,2)==0
        %    Ntime=Ntime;
        % else
        %    Ntime=Ntime+1;
        % end
        % time=delta_t*(1:t_end);
        
        
        %Signal émis
        for ii=1:Ntime
            if time(ii)<=3/Fc
                signal(ii)=1/2*sin(omc*time(ii))*(1-cos(omc*time(ii)/3));
            end
        end
        %cut the trailing zeros in signal
        % signal = signal(1:find(abs(signal) > 0,1,'last')) ;
        %add a few trailing zeros
        signal = [signal zeros(1,100)] ;
        
        Ntime = length(signal) ;
        frequence=(1/delta_t)/Ntime*(0:Ntime-1);
        
        Nfmin=find(frequence<=Fmin);
        Nfmax=find(frequence>=Fmax);
        spectrum=fft(signal);
        Amp = 'f [N]'
        switch(figure_switch)
            case{'display'}%représentation du signal émis et de son spectre
                figure
                subplot(2,1,1)
                plot(time(1:length(signal)),signal)
                subplot(2,1,2)
                plot(frequence,abs(spectrum))
                axis([Fmin*0.2 Fmax*1.8 0 max(abs(fft(signal)))])
        end
    case{'GaussDC', 'gausshape' }
        
        %define the length of the source signal
        
        Ntime=length(time);
        
        signal_length = Source.signal_length  ;%desired length of point source excitation
        steps = signal_length/delta_t    ;%resulting timesteps
        steps = round(steps) ;
        signal_length = (steps)*delta_t   ;%real source length
        
        if isequal(type_switch, 'GaussDC') == 1
            
            %   rise_time = round((Source.Trise)*  steps)  ;
            %   hold_time =  steps - rise_time;
            
            %   if mod(rise_time,2) == 0
            %       rise_time = rise_time +1;
            %      steps = steps + 1 ;
            %      signal_length = (steps)*delta_t  ; %real source length
            %  end
            %  M0_diff_max = Source.Magnitude/floor(rise_time/2);   ; %maximum force to be emitted
            %  delta_force = M0_diff_max/steps   ;
            %  slope1=linspace(0,M0_diff_max,ceil(0.5*rise_time)) ;
            %  slope2 = fliplr(slope1(1:end-1))  ;
            %  GaussDC = [slope1 slope2 zeros(1,hold_time)] ;
            %            signal = cumsum(GaussDC) ;
            x = 0:signal_length/steps:signal_length;
            %             signal = gaussmf(x,[x(floor(0.1*length(x))) x(ceil(0.5*length(x)))]);
            if floor(Source.Trise*length(x)) == 0
                GaussDC = [gaussmf(x,[x(2) x(ceil(0.5*length(x)))])];
                
            else
                GaussDC = [gaussmf(x,[x(floor(Source.Trise*length(x))) x(ceil(0.5*length(x)))])];
            end
            signal = Mo*cumsum(GaussDC).*delta_t ;
            Amp = 'u [m]'
        elseif isequal(type_switch, 'gausshape') == 1
            %             don't start at 0 to avoid problems with
            x = 0:signal_length/steps:signal_length;
            %             signal = gaussmf(x,[x(floor(0.1*length(x))) x(ceil(0.5*length(x)))]);
            if floor(0.1*length(x)) == 0
                signal = [gaussmf(x,[x(2) x(ceil(Source.Trise*length(x)))])];
                
            else
                signal = Mo*[gaussmf(x,[x(floor(0.1*length(x))) x(ceil(Source.Trise*length(x)))])];
            end
            Amp = 'f [N]'
            
        end
        
    case{'ramp'}
        
        %define the length of the source signal
        Ntime=length(time);
        signal_length = Source.signal_length  ;%desired length of point source excitation
        steps = signal_length/delta_t    ;%resulting timesteps
        steps = round(steps) ;
        signal_length = (steps)*delta_t   ;%real source length
        
        
        signal =  Mo*(1 + erf(...
            (0:signal_length/steps:signal_length)/Source.Trise)) ; %moment function
        Amp = 'u [m]'
        
end



signal = Source.Sign.*signal;
switch(type_switch)
    case{'sin_mod','ramp','gausshape' }
        deriv = diff(signal)./delta_t;
        SignalInt = cumsum(signal).*delta_t;
    case{'GaussDC'}
        deriv = diff(signal)./delta_t ;
        SignalInt = cumsum(signal).*delta_t
end

% PlotInt = SignalInt*SigIntFac ;
Plotderiv = deriv ;
PlotInt = SignalInt ;
% define labels depending on force
switch(type_switch)
    case{'sin_mod','gausshape' }
        LabelsSignal = '$\frac{\partial F}{\partial t}$';
        LabelsDeriv =  ['$\frac{d}{d t} F$'];
%         LabelsInt =  ['$(F(t_0) - \int F(t) \,dt)$' ];
        LabelsInt =  ['$F(t_0) - \sum_{t_0}^t F(t) \Delta t$' ];

        
    case{'ramp', 'GaussDC'}
        LabelsSignal = '$u_x$';
        LabelsDeriv =  ['$\frac{d}{d t} u_x$'];
end
L = length(signal);

switch(figure_switch)
    
    case{'display','save'}%
        figure
        
        %         xdft = spectrum(1:NFFT/2+1);
        %         psdx = (1/((1/delta_t)*NFFT)) * abs(xdft).^2;
        %         psdx(2:end-1) = 2*psdx(2:end-1);
        %         dB_PSD = 10*log10(psdx) ;
        
        %         xdft = spectrum2(1:NFFT/2+1);
        %         psdx = (1/((1/delta_t)*NFFT)) * abs(xdft).^2;
        %         psdx(2:end-1) = 2*psdx(2:end-1);
        %         dB_PSD2 = 10*log10(psdx) ;
        %
                
        subplot(2,2,2)
        %         plot(f,(dB_PSD))
        %         hold on;
        %         plot(f,(dB_PSD2))
        %         ylabel('PSD $[\frac{dB}{Hz}]')
        
        if exist('LabelsInt') == 0
            plot(time(1:length(Plotderiv))*1000,Plotderiv)
%             legend( LabelsDeriv)
            
        else
            plot(time(1:length(Plotderiv))*1000,-fliplr(PlotInt(1:end-1)))
%             legend( LabelsInt)
        end
                ylabel('$\frac{\partial M_0}{\partial t} \; [\frac{Nm}{s}]$')

        legend('Derivative')



        grid on
        subplot(2,2,1) ;
        PlotSignal = signal ; %normalised_diff(signal,100)
        plot(time(1:length(PlotSignal))*1000,PlotSignal); hold on
        xlabel('t $[ms]$'); hold on
%         switch(type_switch)
%             case{'sin_mod','gausshape' }          
%                 ylabel('$Amplitude \; [\frac{N}{s}]$')
%             case{'ramp', 'GaussDC'}
%                 ylabel('$Amplitude\; [mm]$');
% %                 LabelsDeriv =  ['$Amplitude \; [\frac{mm}{s}]' '$'];
%         end
%         ylabel('$M_0 \; [Nm]$')
        ylabel('$M_0 \; [Nm]$')

        legend('Physical source','Location','northwest')

        grid on

        
        
        
        
        
        
        xlabel('t $[ms]$'); hold on
%         switch(type_switch)
%             case{'sin_mod','gausshape' }
%                 ylabel('$Amplitude  \; [N]$')
%             case{'ramp', 'GaussDC'}
% %                 LabelsSignal = '$u_x \; [mm]$';
%                 ylabel(['$Amplitude \; [\frac{mm}{s}]' '$']);
%         end
        subplot(2,2,[3 4])  ; hold on
        PosOne = find(Source.sourcesvector == 1) ;
        
        switch(Source.displaytype)
            case{'speed'}
                % time
                % figure; plot([1:length(Source.speedoverdistance)].*delta_t,...
                %         Source.speedoverdistance')
                % space
                speedoverdistance = Source.sourcesvector(PosOne).*...
                    (Params.Trans.spacingMm/1000)./[diff(PosOne)*delta_t diff(PosOne(end-1:end)*delta_t)] ;
                Source.speedoverdistance = zeros(size(Source.sourcesvector));
                for ii = 1:length(Source.sourcesvector)
                      if ismember(ii,PosOne)
                        Source.speedoverdistance(ii) = speedoverdistance(find(PosOne== ii)) ;
                    elseif ii > 1
                        Source.speedoverdistance(ii) = Source.speedoverdistance(ii-1)
                    else
                        Source.speedoverdistance(ii) = Source.speedoverdistance(ii) 
                    end
                end
                Trans = evalin('base','Trans')
    
                Source.spacevector = cumsum((Source.speedoverdistance.*(delta_t))*1000);
                plot(Source.spacevector(PosOne)-( Source.spacevector(PosOne(end)) - Trans.spacingMm*128),...
                    Source.speedoverdistance(PosOne)')
                xlabel('$x \;[mm]$')
                ylabel('$c_r [\frac{m}{s}]$')
%                 ylim([mean(Source.speedoverdistance)-1,max(Source.speedoverdistance)+1])
                legend('rupture speed')
                grid on
            case{'density'}
                %         f = f./1000 ;
                %         set(gca,'XTickLabel',sprintf('%7.0f\n',[round(f(1)) round(f(3)) round(f(10))]  ));
                subplot(2,2,[3 4])  ; hold on
                for ii = 1:length(PosOne)
                    plot((PosOne(ii) : PosOne(ii)+length(signal)-1).*delta_t*1000,...
                        signal,'b')
                end
                xlabel('t [ms]')
                ylabel(LabelsSignal)
        end
        
    case{'save'}
        
        saveas(gcf,Params.filename,'svg')
        %         saveas(gcf,Params.filename,'fig')
%         cleanfigure;matlab2tikz( [Params.filename '.tex'], 'height', '\fheight', 'width', '\fwidth' )
%         
end
%%
% MaxDeriv = max(abs(deriv(:))); MaxSig = max(abs(signal(:)));MaxInt = max(abs(SignalInt(:)))
% SigDerivFac = MaxSig/MaxDeriv ;
% SigIntFac =  MaxSig/MaxInt ;
% Plotderiv = deriv*SigDerivFac ;
% PlotInt = SignalInt*SigIntFac ;
% 
% % define labels depending on force
% switch(type_switch)
%     case{'sin_mod','gausshape' }
%         LabelsSignal = '$F \; [N]$';
%         LabelsDeriv =  ['$\frac{d}{d t} F \; [\frac{N}{s}]*' num2str(round(SigDerivFac,1)) '$'];
%         LabelsInt =  ['$(F(t_0) - \int F(t) \,dt)*' num2str(round(SigIntFac)) '\; [{N}{s}]$'];
%         
%     case{'ramp', 'GaussDC'}
%         LabelsSignal = '$u_x \; [mm]$';
%         LabelsDeriv =  ['$\frac{d}{d t} u_x \; [\frac{mm}{s}]*' num2str(round(SigIntFac,1)) '$'];
% end
% L = length(signal);
% % NFFT = 2^nextpow2(L) ;
% % f = (1/delta_t)/2*linspace(0,1,NFFT/2+1);    %(1:NFFT/2+1)
% [f, spectrum] = spect(signal, 1/delta_t,'Plot',0,'PowerTwo' ,1)
% [f2, spectrum2] = spect(Plotderiv, 1/delta_t,'Plot',0,'PowerTwo' ,1)
% [f3, spectrum3] = spect(PlotInt, 1/delta_t,'Plot',0,'PowerTwo' ,1)
% 
% % spectrum=fft(signal,NFFT)/L;
% % spectrum2=fft(Plotderiv,NFFT)/L;
% % spectrum3=fft(PlotInt,NFFT)/L;
% 
% switch(figure_switch)MaxDeriv = max(abs(deriv(:))); MaxSig = max(abs(signal(:)));MaxInt = max(abs(SignalInt(:)))
% SigDerivFac = MaxSig/MaxDeriv ;
% SigIntFac =  MaxSig/MaxInt ;
% Plotderiv = deriv*SigDerivFac ;
% PlotInt = SignalInt*SigIntFac ;
% 
% % define labels depending on force
% switch(type_switch)
%     case{'sin_mod','gausshape' }
%         LabelsSignal = '$F \; [N]$';
%         LabelsDeriv =  ['$\frac{d}{d t} F \; [\frac{N}{s
%     
%     case{'display','save'}%
%         figure
%         
%         %         xdft = spectrum(1:NFFT/2+1);
%         %         psdx = (1/((1/delta_t)*NFFT)) * abs(xdft).^2;
%         %         psdx(2:end-1) = 2*psdx(2:end-1);
%         %         dB_PSD = 10*log10(psdx) ;
%         
%         %         xdft = spectrum2(1:NFFT/2+1);
%         %         psdx = (1/((1/delta_t)*NFFT)) * abs(xdft).^2;
%         %         psdx(2:end-1) = 2*psdx(2:end-1);
%         %         dB_PSD2 = 10*log10(psdx) ;
%         %
%         
%         subplot(2,2,1) ;
%         PlotSignal = signal ; %normalised_diff(signal,100)
%         plot(time(1:length(PlotSignal))*10000,PlotSignal); hold on
%         xlabel('t $[ms \times 10^{-1}$]'); hold on
%         ylabel('Amplitude')
%         
%         if exist('LabelsInt') == 0
%             plot(time(1:length(Plotderiv))*10000,Plotderiv)
%             legend(LabelsSignal, LabelsDeriv)
%             
%         else
%             plot(time(1:length(Plotderiv))*10000,fliplr(PlotInt(1:end-1)))
%             legend(LabelsSignal, LabelsInt)
%         end
%         
%         subplot(2,2,2)
%         %         plot(f,(dB_PSD))
%         %         hold on;
%         %         plot(f,(dB_PSD2))
%         %         ylabel('PSD $[\frac{dB}{Hz}]')
%         if exist('LabelsInt') == 0
%             
%             plot(f,abs(spectrum)); hold on
%             plot(f2,abs(spectrum2))
%             legend(LabelsSignal, LabelsDeriv)
%             freqcutoff = find(spectrum>=(0.01*(max(abs(spectrum)))),1,'last');
%             
%         else
%             plot(f,abs(spectrum)); hold on
%             plot(f3,abs(spectrum3))
%             legend(LabelsSignal, LabelsInt)
%             freqcutoff = find(spectrum3>=(0.01*(max(abs(spectrum3)))),1,'last');
%             
%         end
%         if isfield(Source,'freqcutoff') == 0
%             Source.freqcutoff =f(freqcutoff);
%         end
%         
%         
%         ylabel('$|P(f)|$')
%         xlabel('f [Hz]')
%         xlim([0 Source.freqcutoff])
%         subplot(2,2,[3 4])  ; hold on
%         PosOne = find(Source.sourcesvector == 1) ;
%         
%         switch(Source.displaytype)
%             case{'speed'}
%                 %% time
%                 % figure; plot([1:length(Source.speedoverdistance)].*delta_t,...
%                 %         Source.speedoverdistance')
%                 %% space
%                 speedoverdistance = Source.sourcesvector(PosOne).*...
%                     (Params.Trans.spacingMm/1000)./[diff(PosOne)*delta_t diff(PosOne(end-1:end)*delta_t)] ;
%                 Source.speedoverdistance = zeros(size(Source.sourcesvector));
%                 for ii = 1:length(Source.sourcesvector)
%                     if ismember(ii,PosOne)
%                         Source.speedoverdistance(ii) = speedoverdistance(find(PosOne== ii)) ;
%                     elseif ii > 1
%                         Source.speedoverdistance(ii) = Source.speedoverdistance(ii-1)
%                     else
%                         Source.speedoverdistance(ii) = Source.speedoverdistance(ii)
%                     end
%                 end
%                 Source.spacevector = cumsum((Source.speedoverdistance.*(delta_t))*1000);
%                 plot(Source.spacevector(PosOne),...
%                     Source.speedoverdistance(PosOne)')
%                 xlabel('$x [cm]$')
%                 ylabel('$c_r [\frac{m}{s}]$')
%                 grid on
%             case{'density'}
%                 %         f = f./1000 ;
%                 %         set(gca,'XTickLabel',sprintf('%7.0f\n',[round(f(1)) round(f(3)) round(f(10))]  ));
%                 subplot(2,2,[3 4])  ; hold on
%                 for ii = 1:length(PosOne)
%                     plot((PosOne(ii) : PosOne(ii)+length(signal)-1).*delta_t*1000,...
%                         signal,'b')
%                 end
%                 xlabel('t [ms]')
%                 ylabel(LabelsSignal)
%             case{'point'}
%                 clf
%                 
%                 %         xdft = spectrum(1:NFFT/2+1);
%                 %         psdx = (1/((1/delta_t)*NFFT)) * abs(xdft).^2;
%                 %         psdx(2:end-1) = 2*psdx(2:end-1);
%                 %         dB_PSD = 10*log10(psdx) ;
%                 
%                 %         xdft = spectrum2(1:NFFT/2+1);
%                 %         psdx = (1/((1/delta_t)*NFFT)) * abs(xdft).^2;
%                 %         psdx(2:end-1) = 2*psdx(2:end-1);
%                 %         dB_PSD2 = 10*log10(psdx) ;
%                 %
%                 
%                 subplot(2,1,1) ;
%                 PlotSignal = signal ; %normalised_diff(signal,100)
%                 plot(time(1:length(PlotSignal))*1000,PlotSignal); hold on
%                 xlabel('t $[ms]$'); hold on
%                 switch(type_switch)
%                     case{'sin_mod','gausshape' }
%                         ylabel('$Amplitude \; [\frac{N}{s}]$')
%                     case{'ramp', 'GaussDC'}
%                         ylabel('$Amplitude\; [mm]$');
%                         %                 LabelsDeriv =  ['$Amplitude \; [\frac{mm}{s}]' '$'];
%                 end
%                 legend(LabelsSignal)
%                 grid on
%                 
%                 subplot(2,1,2)
%                 %         plot(f,(dB_PSD))
%                 %         hold on;
%                 %         plot(f,(dB_PSD2))
%                 %         ylabel('PSD $[\frac{dB}{Hz}]')
%                 
%                 if exist('LabelsInt') == 0
%                     plot(time(1:length(Plotderiv))*1000,Plotderiv)
%                     legend( LabelsDeriv)
%                     
%                 else
%                     plot(time(1:length(Plotderiv))*1000,fliplr(PlotInt(1:end-1)))
%                     legend( LabelsInt)
%                 end
%                 grid on
%                 
%                 
%                 figure
%                 %         xdft = spectrum(1:NFFT/2+1);
%                 %         psdx = (1/((1/delta_t)*NFFT)) * abs(xdft).^2;
%                 %         psdx(2:end-1) = 2*psdx(2:end-1);
%                 %         dB_PSD = 10*log10(psdx) ;
%                 %         xdft = spectrum2(1:NFFT/2+1);
%                 %         psdx = (1/((1/delta_t)*NFFT)) * abs(xdft).^2;
%                 %         psdx(2:end-1) = 2*psdx(2:end-1);
%                 %         dB_PSD2 = 10*log10(psdx) ;
%                 %
%                 subplot(1,2,1) ;
%                 PlotSignal = signal ; %normalised_diff(signal,100)
%                 plot(time(1:length(PlotSignal))*1000,PlotSignal); hold on
%                 xlabel('t $[ms]$'); hold on
%                 switch(type_switch)
%                     case{'sin_mod','gausshape' }
%                         ylabel('$Amplitude \; [\frac{N}{s}]$')
%                     case{'ramp', 'GaussDC'}
%                         ylabel('$Amplitude\; [mm]$');
%                         %                 LabelsDeriv =  ['$Amplitude \; [\frac{mm}{s}]' '$'];
%                 end
%                 legend(LabelsSignal)
%                 grid on
%                 subplot(1,2,2)
%                 %         plot(f,(dB_PSD))
%                 %         hold on;
%                 %         plot(f,(dB_PSD2))
%                 %         ylabel('PSD $[\frac{dB}{Hz}]')
%                 if exist('LabelsInt') == 0
%                     plot(time(1:length(Plotderiv))*1000,Plotderiv)
%                     legend( LabelsDeriv)
%                 else
%                     plot(time(1:length(Plotderiv))*1000,fliplr(PlotInt(1:end-1)))
%                     legend( LabelsInt)
%                 end
%                 grid on
%         end
%         
%     case{'save'}
%         
%         saveas(gcf,Params.filename,'svg')
%         %         saveas(gcf,Params.filename,'fig')
%         cleanfigure;matlab2tikz( [Params.filename '.tex'], 'height', '\fheight', 'width', '\fwidth' )
%         
% end

%%
end
