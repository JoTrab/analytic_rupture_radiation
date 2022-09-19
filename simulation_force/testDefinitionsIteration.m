clear all

list1 = {}
list1{1} = 'Vp'
list1{2} = 'Dpl'
for ii = 1:2
clear F_all_t Ut vt Fieldterms
 Params.Vp_Dpl = list1{ii}
addpath(genpath('~/work/functions/my_functions/Matlab/elastic_wavefield_analytical/simulation_force/'))
   % decide if you want displacement or Vp as primary output
Trans.numelements_true = 128;
Trans.size_factor = 10; %reduce size and increase distance in z
Trans.Depthpoints = 1500/Trans.size_factor;  %delta depth is 1/8 of ...
        % wavelength ultrasound scaled by size factor
Trans.frequency = 5.208; % nominal frequency in MHz
Trans.Bandwidth = [4, 7] ;
Trans.elementWidth = .250; % width in mms
Trans.spacingMm = .298;   % Spacing between elements in mm.
Params.Fs = 3000;
% Params.returnSwitch = 'Field' ; % 'Green', 'Field' ,'none'
Params.Field_only = 'off'

Rupture.rup_speed_start_factor = 2  %1.0 ;
Rupture.rup_speed_end_factor = 2    % 0.1  ;%sqrt(2)  % sqrt(2) ;  %sqrt(2)
%Give a maximum allowed rupture time
Rupture.max_rup_time = 0.03  ;  %total time from one end of rupture to the other
%Give a minimum allowed rupture time
Rupture.min_rup_time = 0.000001  ;
% Rupture.max_rup_steps = 250; %maximum time_points wished for calculation
Rupture.length_total = 0.003  ;% in m accelerating and steady length
%add a constant speed vector at start/end
Rupture.final_speed = 0.003 ; %add number of sources / points/ timesteps in meter that work at final speed
Rupture.start_speed = 0.0000001

Source.increase_length = 20   ;%5 ; %increases number of  timesteps since it
    %decreases delta_t to allow for stretching
    %becomes important if the source length becomes too small due too many
    %shifted sources. d_shift min defines minimal delta_t before scalintg by this factor
Params.observation_dir = 'vertical' ;   %'vertical' 'horizontal'  Gx - horizontal  Gy - side Gz - vertical
                                                % vertical is z and horizontal is x  , side is y    % defines the
                                                         % position of the probe and thus which vp is retrievec, vx,vy,vz

                                                         
                                                         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%SOURCE DIRECTION OBSERVATION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                      
Source.Direction = 90  
Rupture.rup_direction = 'pos' ; %'pos' 'neg'  pos is right or up, neg is down or left

%define relative to observation direction 90 means right for vertical ,
%  0 means right for horizontal; define always positive
                                                   
Params.observation_plane = 'x-z' % 'x-y'  ,'x-z' %defines the plane of observation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Params.s_speed = 6.9;
Params.p_speed = 1540;
Params.density = 1000;
Params.t_end = 0.01    ;% 0.025 ; %final simulation time  in seconds  - careful that not too short!!
Params.calc_type = 'discrete'  ;% % smooth should be obsolete!! 'smooth' 'discrete'
Params.video = 'on'  ;  %'on'  %%'off'
                                %Source point also automatically defines the position
                                %2 source points: above or on the side top to bottom x-z plane
                                %3 source points: x-y plane

                                %%Source.Point = [floor(0.5*Trans.Depthpoints) floor(0.3*Trans.numelements_true) 0.003];

%%%%%%%%%%%%%% SOURCE POINT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Source.Point = [Trans.Depthpoints floor(0.3*Trans.numelements_true)];
                                    %---> source length is calculated from this.
                                    %only use positive non-values of the source as source vector
                                    %time sampling equal to rupture time of 1 point to the next:
                                    %Define a more or less centred SourcePoint
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Source.Sign = 1   ; %Direction of the Source.
% Source.Cycles = 1;
% Source.Frequency_Hz = 250;
Source.Magnitude = 1;  %M0 for the double couple
% Source.Bandwidth = 0.6; % bandwidth of gausspuls
%for the case of a non-accelerating source ; for accelerating source: adapt
Source.type_switch = 'gausshape'   ;
Source.figure_switch = 'save'  ;  %'display', 'no' ,'save'
Source.Trise = 0.5 ;  %defines tris in ramp. make it 0.5 for gaussian
Source.d_shift_min = 0.5*Source.Trise ; % start  or end depends if accel or decel
% Source.d_shift_start = 2*Source.Trise ; %

Params.filename = strcat('Force_', num2str(Rupture.rup_speed_start_factor), '_', ... 
    num2str(Rupture.rup_speed_end_factor), '_', Params.observation_dir, '_',  ... 
    Params.observation_plane, '_', ...
    num2str(Params.calc_type), '_', Params.Vp_Dpl )  ;

%Params.filename  = 'test'
% obsolete ? Params.Field_only = 'off' %stop before shifting the sources
Params.testrun = 'off'  ;  %'on', 'off'





%%
if isequal(Params.testrun,'on') == 1
    [~,~,~,Source, Params, Rupture] =  wrapper_function(Trans, Source, Params,  Rupture) ;
else
    if isequal(Params.Field_only,'on') 
        [ ~,Fieldterms,~,Source, Params, Rupture] = wrapper_function(...
        Trans, Source, Params,  Rupture)   ;
    else
        [Ut, Fieldterms, vt, Source, Params, Rupture,F_all_t] = wrapper_function(...
        Trans, Source, Params,  Rupture)   ;
    end
end
end