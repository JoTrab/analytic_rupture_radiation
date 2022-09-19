
%% Define the space time ( in points  ) vector
%the maximum time distance between to rupture points (the beginning)
%get the length of the acceleration vector without adapting delta_t
clear PosOne
if rup_speed_end > rup_speed_start
    rup_speed_max = rup_speed_end ;
    rup_speed_min = rup_speed_start ;
else
   rup_speed_max = rup_speed_start ;
   rup_speed_min = rup_speed_end ;
end
timestepping_factor_max =  round(rup_speed_max / rup_speed_min);
source_vector_temp = zeros(1,sum(1:(timestepping_factor_max-1))) ;  %+timestepping_factor_max+2)

%% create a vector in time for the source points

for i = 1: (timestepping_factor_max )  % take the amount by which delta_t is smaller into account
    if i == 1
        PosOne(i) = 1 ; %set the first point to 1 --> rupture
    else
        PosOne(i)= PosOne(i-1) + timestepping_factor_max-(i-1)  ;
    end
end

% first_length = length(PosOne)
%adapt delta_t using increase_acc_length
increase_acc_length = Rupture.accelerating_length_points - length(PosOne) ;
if  increase_acc_length > 0
    increase_length_factor = (  round(increase_acc_length / (rup_speed_max / rup_speed_min)) -1  ) +1   ;
    clear PosOne
    %add the user chosen increase of source_length
    % increase_length_factor = increase_length_factor * Source.increase_length;
    timestepping_factor_max = round((increase_length_factor) * (rup_speed_max / rup_speed_min)); %distance between ones in the
%     source_vector_temp = zeros(1,sum(1:(timestepping_factor_max-1))) ;  %+timestepping_factor_max+2)
    for i = 1: (timestepping_factor_max - (increase_length_factor-1))  % take the amount by which delta_t is smaller into account
        if i == 1
            PosOne(i) = 1 ; %set the first point to 1 --> rupture
        else
            PosOne(i)= PosOne(i-1) + timestepping_factor_max-(i-1)  ;
        end
    end
else
    increase_length_factor = 1 ;
%     source_vector_temp =zeros(1,PosOne(end)) ;
end

source_vector_temp =zeros(1,PosOne(end)) ;
source_vector_temp(PosOne) = 1  ;
if rup_speed_end < rup_speed_start
  %flip vor decelarating
  source_vector_temp = fliplr(source_vector_temp);
  time_diff_end = timestepping_factor_max ;
  time_diff_start = increase_length_factor ;
else
  time_diff_start = timestepping_factor_max ;
  time_diff_end = increase_length_factor ;
end


%correct acelerating length
Rupture.accelerating_length_points = length(PosOne);
% correct true add_rupEnd numer
% add_rupBoth_source_nr = Rupture.rup_length_points -  Rupture.accelerating_length_points ;
% add_rupStart_source_nr = ceil(Rupture.add_rupStart_source_nr/Rupture.add_rupStart_source_nr * add_rupBoth_source_nr)
% add_rupEnd_source_nr = add_rupBoth_source_nr - add_rupStart_source_nr

% Rupture.add_rupEnd_source_nr = add_rupEnd_source_nr ;
% Rupture.add_rupStart_source_nr = add_rupStart_source_nr;
add_rupEnd_source_nr = Rupture.add_rupEnd_source_nr ;
add_rupStart_source_nr = Rupture.add_rupStart_source_nr ;


add_rupEnd = zeros(add_rupEnd_source_nr*time_diff_end,1)'; %adds a certain nu,ber of sources at end speed
add_rupStart = zeros(add_rupStart_source_nr*time_diff_start,1)'; %adds a certain nu,ber of sources at end speed

add_rupEnd(time_diff_end:time_diff_end:end) =1
add_rupStart(time_diff_start:time_diff_start:end) =1
add_rupStart = fliplr(add_rupStart) ;

lastOne = find(source_vector_temp,1,'last') ;
source_vector_temp = [ add_rupStart source_vector_temp(1:lastOne) add_rupEnd] ;
nr_sources = length(PosOne) + add_rupEnd_source_nr + add_rupStart_source_nr;%source space length
%redefine PosOne after adding the coSnstant speed values
PosOne = find(source_vector_temp == 1);


%stretch the vector for smoother calculation ?
PosOne(1:end) = PosOne(1:end) * Source.increase_length ;
source_vector_temp = zeros(1,sum(1:(timestepping_factor_max-1)) * Source.increase_length ) ;
source_vector_temp(PosOne) = 1 ;
rup_steps = length(source_vector_temp) ; %source time  length

%% define time , delta time and SOurce time and length (in seconds)

%delta_t =   (1/increase_length_factor) * ((Trans.spacingMm/1000) / rup_speed_end) ;
% The minimum shift time in seconds
if length(PosOne)>1
    Source.d_shift_min_points =  min(PosOne(end)- PosOne(end-1), PosOne(2) - PosOne(1)) ;
else
    Source.d_shift_min_points = 0 ;
    % in the case of one source point return to wrapper
    display('One Source point only')
end

%caulation needs to be done on bigger grid (due to shift)
if Source.d_shift_min_points ~= 0 
    Trans.numelements = Trans.numelements_true +Trans.numelements_true; % needed to avoid wraparound and keep source at defined position
else
    Trans.numelements = Trans.numelements_true; % needed to avoid wraparound and keep source at defined position
end

Trans.ElementPos = zeros(Trans.numelements,4);
Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));


if Source.d_shift_min_points ~= 0
    % % delta_t = Rupture.total_rup_time/PosOne(end) ; %length(source_vector_temp)
    %minimum time for advancing one rupture point
    Source.d_shift_min_seconds = (Trans.spacingMm / 1000) / (max(Rupture.rup_speed_end_factor,Rupture.rup_speed_start_factor) *Params.s_speed);
    delta_t = Source.d_shift_min_seconds/Source.d_shift_min_points ;
    Rupture.total_rup_time = PosOne(end) *delta_t;
    if Rupture.total_rup_time > Rupture.max_rup_time
        %total time from one end of rupture to the other
        error(['rupture time too long by ' num2str(Rupture.total_rup_time-Rupture.max_rup_time) 's --> adapt length or speed'])
    elseif  Rupture.total_rup_time < Rupture.min_rup_time
            error(['rupture time too short by ' num2str(Rupture.min_rup_time-Rupture.total_rup_time) 's --> adapt length or speed'])
    end
    % if isequal(Source.type_switch , 'triangle') == 1
    Source.signal_length = (Source.d_shift_min_seconds) / (Source.d_shift_min);
    % end
else
    reply = input('Rupture length in ms','s');
    if isempty(reply)
      reply = 1;
    end
  
    reply2 = input('timesteps for rupture nucleation','s');
    if isempty(reply2)
      reply2 = 20;
    end  
    
    Source.signal_length = (str2num(reply)/1000)
    delta_t = Source.signal_length/str2num(reply2);
    
    if delta_t > (1/Params.Fs) 
        delta_t = (1/Params.Fs) ;
    end
%     Source.signal_length= 1/Params.Fs ; 
    Rupture.total_rup_time  = Source.signal_length;
end
time= 0:delta_t:t_end;
if length(time) < (Source.signal_length/delta_t)
    error('reduce source signal length')
end

% time = [-fliplr(time(2:10)) time];  %add negative times
total_timesteps = length(time);
timesteps = total_timesteps ;  %timesteps to give in calculation equal to or lower to total timesteps
real_speed_min = (Trans.spacingMm/1000) / (timestepping_factor_max*delta_t)  ;

%linear accelerate on this length ----> scale delta t by increase length
%factor to get the right length












%%

%%define grid

% grid = zeros(Depthpoints,Trans.numelements);
depth_meter = (0.*(1:Depthpoints))' ;
% put the origin in upperleft corner
x_vec = ((Trans.ElementPos(:,1)-Trans.ElementPos(1,1))  ./ 1000)';
%delta_x = x_vec(2) - x_vec(1) ;
x_mat = repmat(x_vec,[length(depth_meter), 1]);
delta_depth = (1./8).*size_factor.*(1540/(Trans.frequency*1e6)) ;

%put Source at bottom
if Source.Point(1) == Depthpoints
    depth_meter(:) = (1:length(depth_meter)).* delta_depth;  %    WRONG!!! due to stelus speckle tracking?
    depth_meter = flipud(depth_meter);
else
    %put Source in medium
    depth_meter(1:Source.Point(1)-1) = fliplr((1: Source.Point(1) -1))  ...
        .* delta_depth;
    depth_meter(Source.Point(1)+1:end) = (1:length(depth_meter) - Source.Point(1)) ...
        .* delta_depth.*(-1);
end
depth_mat = repmat(depth_meter,[1 length(x_vec)]);


if Source.d_shift_min_points ~= 0 
% wave_advance_per_time = (beta*delta_t/(Trans.spacingMm/1000));
    if isequal(Rupture.rup_direction, 'pos') ==1
        source_ind = [Source.Point(1,1), Source.Point(1,2)+Trans.numelements_true]; %shift the source to the right
        %since we need the field calculated on the legt for the shift
    elseif isequal(Rupture.rup_direction, 'neg') ==1
        source_ind = [Source.Point(1,1), Source.Point(1,2)];
        %don t shift the source, the field is already calculated on the
        %right for the shift
    end
else
	if isequal(Rupture.rup_direction, 'pos') ==1
        source_ind = [Source.Point(1,1), Source.Point(1,2)]; %shift the source to the right
        %since we need the field calculated on the legt for the shift
    elseif isequal(Rupture.rup_direction, 'neg') ==1
        source_ind = [Source.Point(1,1), Source.Point(1,2)];
        %don t shift the source, the field is already calculated on the
        %right for the shift
    end
end
Source.Point_input = Source.Point ;
if length(Source.Point_input) ==3
    source_ind = [source_ind Source.Point_input(3)]
end
Source.Point = source_ind;

padnumber_space = nr_sources ; %wave_advance_per_time*nr_sources + 5;
padnumber_time = length(source_vector_temp) + 10 ;%ceil((1/rup_speed_end_factor)*total_timesteps) +5 ;
%padnumber_time = (total_timesteps) +5 ;
padnumber.time = padnumber_time ;
padnumber.space = padnumber_space ;




if isequal(Rupture.rup_direction,'pos') ==1  & Source.Point_input(2)+nr_sources > Trans.numelements_true
    depass = Trans.numelements_true - (Source.Point_input(2)+nr_sources);
    error(['Source exceeds grid by ' num2str(depass) ' points'])
elseif isequal(Rupture.rup_direction,'neg') ==1 & Source.Point_input(2)-nr_sources <= 0
    depass = (Source.Point_input(2)-nr_sources);
    error(['Source exceeds grid by ' num2str(depass) ' points'])
end


% tc = gauspuls('cutoff',Source.Frequency_Hz,Source.Bandwidth,[],-40);  %last truncation of db , 60 % bandwidth
% t = -1*tc : delta_t : 1*tc;  %sample the pulse at 3 khz

%add zeros until the end of calculation to the source term
% source_ux = [source_ux zeros(1, timesteps - length(source_ux))]   ;
% X0 = single(zeros(Depthpoints,Trans.numelements,1));
% X0 = repmat(X0,1,1,timesteps);
% X0 = (zeros(Depthpoints,Trans.numelements,timesteps));
% X1 = X0;
r = zeros(length(depth_meter), length(x_vec));


Params.delta_t = delta_t;
Params.rup_speed = [rup_speed_start rup_speed_end];
Params.beta = beta;
Params.p_speed = p_speed;
Params.time = time;
Params.padnumber = padnumber;
Params.observation_dir = observation_dir;
Params.x_vec = x_vec ;
Params.depth_meter = depth_meter ;
Params.Trans = Trans ;


% Source.X0  = X0              ;
% Source.X1 = X1 ;
Source.Number =  nr_sources ;
Source.sourcesvector = source_vector_temp ;

if isequal(Params.observation_plane,'x-z')
    r(:,:)=sqrt( ( x_mat(:,:) - x_mat(Source.Point(1,1),Source.Point(1,2) )).^2 ...
        + ( abs(depth_mat(:,:) - depth_mat(Source.Point(1,1),Source.Point(1,2) ))).^2) ;
elseif isequal(Params.observation_plane,'x-y')
    r(:,:)=sqrt( ( x_mat(:,:) - x_mat(Source.Point(1,1),Source.Point(1,2) )).^2 ...
    + ( abs(depth_mat(:,:) - depth_mat(Source.Point(1,1),Source.Point(1,2) ))).^2 ...
    +  Source.Point(1,3).^2)  ;
else
    error('Observation plane not defined, define observation plane as x-z or x-y' )
end
% r is here just the distance in the 2d plane. even for the 3d case
clear X0 X1
% [1 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 ]
