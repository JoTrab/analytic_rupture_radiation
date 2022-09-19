Params.observation_dir =   'horizontal' ;   %'vertical' 'horizontal'
Params.s_speed = 6;
Params.p_speed = 1540;
Params.density = 1000;

%% DEfinitions
% check that the convoluted source signal in the near field has a plateau
% ---> source signal long enough
clear all
Trans.numelements_true = 240;
Trans.size_factor = 10; %reduce size and increase distance in z
Trans.Depthpoints = 1000/Trans.size_factor;
Trans.frequency = 5.208; % nominal frequency in MHz
Trans.Bandwidth = [4, 7] ;
Trans.elementWidth = .250; % width in mm
Trans.spacingMm = .298/2;   % Spacing between elements in mm.
Params.Fs = 8000;


Params.Fs = 3000;

Rupture.rup_speed_start_factor = 0.7 ;
Rupture.rup_speed_end_factor = 0.7  ;%sqrt(2)  % sqrt(2) ;  %sqrt(2)
Rupture.max_rup_time = 0.004  ;  %total time from one end of rupture to the other
Rupture.min_rup_time = 0.001  ;
Rupture.max_rup_steps = 250; %maximum time_points wished for calculation
Rupture.length_total = 0.012  ;% in m accelerating and steady length
Rupture.final_speed =0.012 ; %add number of sources / points/ timesteps in meter that work at final speed
Source.increase_length = 20  ;%5 ; %increases number of  timesteps since it decreases delta_t to allow for stretching
%becomes important if the source length becomes too small due too many
%shifted sources. d_shift min defines minimal delta_t before scalintg by this factor
Params.observation_dir =   'horizontal' ;   %'vertical' 'horizontal'
Params.s_speed = 6;
Params.p_speed = 1540;
Params.density = 1000;
Params.t_end = 0.01    ;% 0.025 ; %final simulation time  in seconds  - careful that not too short!!



Params.calc_type = 'discrete'  ;% 'smooth' 'discrete'
Rupture.rup_direction = 'pos' ; %'pos' 'neg'  pos is right or up, neg is down or left
Params.video = 'on'  ;  %'on'  %%'off'
Source.Point = [floor(0.5*Trans.Depthpoints) floor(0.3*Trans.numelements_true)];

%---> source length is calculated from this.


%only use positive non-values of the source as source vector

%time sampling equal to rupture time of 1 point to the next:

%Define a more or less centred SourcePoint
Source.Sign = 1   ; %Direction of the Source.
Source.Direction = 0  %define relative to observation direction 90 means right for vertical , 0 means right for horizontal; define always positive


% Source.Cycles = 1;
% Source.Frequency_Hz = 250;
Source.Magnitude = 1;  %M0 for the double couple
% Source.Bandwidth = 0.6; % bandwidth of gausspuls

%for the case of a non-accelerating source ; for accelerating source: adapt
Source.type_switch = 'gausshape'   ;
Source.figure_switch = 'display'  ;  %'display', 'no'
Source.Trise = 0.5 ;  %defines tris in ramp. make it 0.5 for gaussian
Source.d_shift_min = 0.5*Source.Trise ; % start  or end depends if accel or decel
% Source.d_shift_start = 2*Source.Trise ; %
Params.Field_only = 'on'
addpath(genpath('~/functions/simulation/force/'))
Rupture.rup_direction = 'pos' ; %'pos' 'neg'  pos is right or up, neg is down or left

%%  RIGHT
Source.Point = [floor(0.5*Trans.Depthpoints)-1 floor(0.3*Trans.numelements_true)];
Source.Sign = 1   ; %Direction of the Source.
Source.Direction = 0
Params.filename = [Params.observation_dir Rupture.rup_direction num2str(Source.Direction) num2str(Source.Sign)  ] ;
Params.testrun = 'off'  ;
if isequal(Params.testrun,'on') == 1
    [~,~,~,Source, Params, Rupture] =  wrapper_function(Trans, Source, Params,  Rupture) ;
else
    [~, Fieldterms, ~, Source, Params, Rupture ] = wrapper_function(Trans, Source, Params,  Rupture)   ;

end
Fieldterms_right = Fieldterms.Field_all;
clear Source.X0 Source.X1  Source.X2

%% LEFT
Source.Point = [floor(0.5*Trans.Depthpoints)+1 floor(0.3*Trans.numelements_true)];

Source.Sign = -1   ; %Direction of the Source.
Source.Direction = 0  %define relative to observation direction 90 means right for vertical , 0 means right for horizontal; define always positive
Params.filename = [Params.observation_dir Rupture.rup_direction num2str(Source.Direction) num2str(Source.Sign)  ] ;

if isequal(Params.testrun,'on') == 1
    [~,~,~,Source, Params, Rupture] =  wrapper_function(Trans, Source, Params,  Rupture) ;
else
    [~, Fieldterms, ~, Source, Params, Rupture ] = wrapper_function(Trans, Source, Params,  Rupture)   ;

end
Fieldterms_left = Fieldterms.Field_all;
clear Source.X0 Source.X1  Source.X2

%% DOWN
Source.Point = [floor(0.5*Trans.Depthpoints) floor(0.3*Trans.numelements_true)-1];
Source.Sign = -1   ; %Direction of the Source.
Source.Direction = 90
Params.filename = [Params.observation_dir Rupture.rup_direction num2str(Source.Direction) num2str(Source.Sign)  ] ;

if isequal(Params.testrun,'on') == 1
    [~,~,~,Source, Params, Rupture] =  wrapper_function(Trans, Source, Params,  Rupture) ;
else
    [~, Fieldterms, ~, Source, Params, Rupture ] = wrapper_function(Trans, Source, Params,  Rupture)   ;

end
Fieldterms_down = Fieldterms.Field_all;
clear Source.X0 Source.X1  Source.X2

%% UP
Source.Point = [floor(0.5*Trans.Depthpoints) floor(0.3*Trans.numelements_true)+1];
Source.Sign = 1   ; %Direction of the Source.
Source.Direction = 90
Params.filename = [Params.observation_dir Rupture.rup_direction num2str(Source.Direction) num2str(Source.Sign)  ] ;

if isequal(Params.testrun,'on') == 1
    [~,~,~,Source, Params, Rupture] =  wrapper_function(Trans, Source, Params,  Rupture) ;
else
    [~, Fieldterms, ~, Source, Params, Rupture ] = wrapper_function(Trans, Source, Params,  Rupture)   ;

end
clear Source.X0 Source.X1  Source.X2
Fieldterms_up = Fieldterms.Field_all;



%% shift sources
Field = Fieldterms_up + Fieldterms_down + Fieldterms_left + Fieldterms_right ;

% Field_2 = Field(:,:,1:floor(1/Params.Fs  /Params.delta_t):size(Field,3)); %low temporal reolution
% vt_Field2 = diff(Field_2(:,:,:),[],3)./ (1/Params.Fs)  ;
% vt_Field = diff(Field_high(:,:,:),[],3)./delta_t ;%./ (1/Params.Fs)  ;
'finish Field'
%%

rup_steps = Params.rup_steps;
padnumber_space = Params.padnumber.space;
padnumber_time = Params.padnumber.time;
PosOne = Params.PosOne;
tic
shift_sources %gives Field_pad
toc
'finish shift'

%% particle velocity

% vt(:,:,:) = diff(Field(:,:,:),[],3) ./ Params.delta_t   ;
% subplot(3,1,3);imagesc(squeeze(Field_static(:,:)),limits);colorbar;title(['Gstat' num2str(i)]);

%vt(:,:,:) = diff(Field_pad(:,:,:),[],3) ./ Params.delta_t   ;
F_all_t =Field_pad;
Field_pad = Field_pad(:,:,1:floor(1/Params.Fs  /Params.delta_t):size(Field_pad,3)-Params.padnumber.time);
% Field = Field(:,:,1:floor(1/Params.Fs  /Params.delta_t):size(Field,3));
% % field is cleared
vt(:,:,:) = diff(Field_pad(:,:,:),[],3) ./ (1/Params.Fs)  ;

'finish vt'

%% save Data
clear Source.Field_padded
%cd /media/johnny/sata_data/Simulations/vertical/subshear
%save simulation G_fct G_fct_near G_fct_s Field_pad vt Source Params r Trans -v7.3
save([Params.filename '.mat'] , 'Field_pad', 'Fieldterms', 'vt', 'Source', 'Params', 'r', 'Trans', 'Fieldterms_up', 'Fieldterms_down', 'Fieldterms_left', 'Fieldterms_right', '-v7.3')
%% Make VIdeo
limitsdpl = [-0.04*max(max(max(Field_pad))) 0.04*max(max(max(Field_pad)))];
limitsvt  = [-0.04*max(max(max(vt))) 0.04*max(max(max(vt)))];
videoname =  ['DC_' num2str(Rupture.rup_speed_start_factor) '_' num2str(Rupture.rup_speed_end_factor) Params.observation_dir num2str(Rupture.add_rupEnd_source_nr) ]  ;   %Params.filename ;
subplot_pos = '1,2,'
% if isequal(Params.video ,'on') == 1
    Particle_Deplacement_film(Field_pad,vt, Params, [] ,'four_forces', limitsdpl, limitsvt, subplot_pos )
% end
