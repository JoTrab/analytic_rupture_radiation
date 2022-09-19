

% Fc=Source.Frequency_Hz ;
% Fmax=Fc + ceil(0.5*Source.Frequency_Hz*Source.Bandwidth);
% Fmin=Fc - ceil(0.5*Source.Frequency_Hz*Source.Bandwidth);
% omc=2*pi*Fc;

% Tmax=0.1;

Ntemps=length(time);
% if rem(Ntemps,2)==0
%    Ntemps=Ntemps;
% else
%    Ntemps=Ntemps+1;
% end
temps = time;
% temps=delta_t*(1:t_end);
source_signal_length = 0.002  ;%desired length of point source excitation
steps = source_signal_length/delta_t    ;%resulting timesteps
steps = round(steps) ;
source_signal_length = (steps)*delta_t   ;%real source length


max_force = 1   ; %maximum force to be emitted
delta_force = max_force/steps   ;
hold_time =   round((3/4)*  steps)  ; 
rise_time = steps - hold_time; 
if mod(rise_time,2) == 0
    rise_time = rise_time +1;
    steps = steps + 1 ;
    source_signal_length = (steps)*delta_t  ; %real source length
end

    

slope1=linspace(0,max_force,ceil(0.5*rise_time)) ;
slope2 = fliplr(slope1(1:end-1))  ;
triangle = [slope1 slope2 zeros(1,hold_time)] ;
signal = cumsum(triangle) ;



spectrum=fft(signal);
figure
subplot(2,1,1) ;
plot(temps(1:length(signal)),signal); hold on
plot(temps(1:length(triangle)),triangle); 
ylim([0 80]) ;
subplot(2,1,2)
plot(abs(spectrum))
