
function [Ut, Fieldterms , vt, Source, Params, Rupture, varargout] = wrapper_function(Trans, Source, Params,  Rupture)
%%
size_factor = Trans.size_factor ; %reduce size and increase distance in z
Depthpoints = Trans.Depthpoints ;
observation_dir = Params.observation_dir ;
observation_dir = Params.observation_dir ;
p_speed = Params.p_speed ;
density = Params.density;
t_end = Params.t_end ; %final simulation time  in seconds  - careful that not too short!!
calc_type = Params.calc_type ; % 'smooth' 'discrete'
beta = Params.s_speed;

%% define US basics


rup_speed_start_factor = Rupture.rup_speed_start_factor ;
rup_speed_end_factor = Rupture.rup_speed_end_factor ; %sqrt(2)  % sqrt(2) ;  %sqrt(2)
%max_rup_steps = Rupture.max_rup_steps ; %maximum time_points wished for calculation
% add_rupEnd_source_nr = Rupture.add_rupEnd_source_nr  ; %add number of sources / points/ timesteps that work at supershear speed
rup_length_points = round(Rupture.length_total / ((Trans.spacingMm/1000) ));
Rupture.rup_length_points = rup_length_points;
Rupture.add_rupEnd_source_nr = round(Rupture.final_speed  / ((Trans.spacingMm/1000) )); %add number of sources / points/ timesteps in meter that work at final speed
Rupture.add_rupStart_source_nr = round(Rupture.start_speed  / ((Trans.spacingMm/1000) )); %add number of sources / points/ timesteps in meter that work at final speed
Rupture.accelerating_length_points = rup_length_points - ...
    Rupture.add_rupStart_source_nr   - Rupture.add_rupEnd_source_nr;

%%
%%
rup_speed_start = rup_speed_start_factor*beta%  sqrt(2)*beta ;  %0.1*beta ;
rup_speed_end = rup_speed_end_factor * beta ;
%% sim_cluster_define
accelerating_source

% the time point shift calculated in accelerating_source and the Source.d_shift
%define the final Source time length to ensure the shift is done at the right point in the source function



[Source.Function, Source.Function_deriv] = source_signal(Source,Params) ;
% triangle_source_signal
%%
% Source.Function = [zeros(1,10) ones(1,100) zeros(1,10)]
%arrival times s and p
if isequal(Params.testrun,'on') == 1
    Ut = []  ;
    Fieldterms = []  ;
    vt = []  ;
    return
end

Source.X0 = Delta_source_function(time(1:timesteps), r, beta,Source.Point,depth_mat,calc_type)  ;
Source.X1 = Delta_source_function(time(1:timesteps), r, p_speed,Source.Point,depth_mat,calc_type)  ;
Source.X2 = Delta_near_field(r,time,beta,p_speed) ;

wave_type = 's-wave' ;
Green_type_switch ;
wave_type = 'p-wave' ;
Green_type_switch ;
wave_type = 'near-field' ;
Green_type_switch ;
G_fct = G_fct_near + G_fct_s + G_fct_p ;
% Source = rmfield(Source,{'X0','X1','X2'});

Source = rmfield(Source, {'X0','X1', 'X2'}) ;
% clear Source.X0 Source.X1 Source.X2 %G_fct_s G_fct_p G_fct_near
if isequal(Params.Vp_Dpl,'Vp')   % decide if you want displacement or Vp as primary output
    tic
    Field = block_convolution(Source.Function_deriv, G_fct, 60);
    toc
elseif isequal(Params.Vp_Dpl,'Dpl')
    tic
    Field = block_convolution(Source.Function, G_fct, 60);
    toc
end
Params.rup_steps = rup_steps;
Params.PosOne = PosOne;
Fieldterms.Field_all = Field;  %high temporal reolution
Fieldterms.G_all = G_fct;
Fieldterms.G_fct_s = G_fct_s;
Fieldterms.G_fct_near = G_fct_near;


% % % % Field = Fieldterms.G_fct_s ;

if isequal(Params.Field_only, 'on') == 1
    Ut = [];
    vt = [];
    F_all_t = [];
    RadPatFigure
    
    return
end
% Field_2 = Field(:,:,1:floor(1/Params.Fs  /Params.delta_t):size(Field,3)); %low temporal reolution
% vt_Field2 = diff(Field_2(:,:,:),[],3)./ (1/Params.Fs)  ;
% vt_Field = diff(Field_high(:,:,:),[],3)./delta_t ;%./ (1/Params.Fs)  ;
'finish Field'
%%
% no padding if only one source

% delete green functions and fieldterms

clear G_fct_s G_fct_near G_fct


if length(find(Source.sourcesvector==1)) > 1
    tic
    shift_sources %gives Ut
    toc
    'finish shift'
else
    Field_pad = Field ;
end


if Source.ReleaseTogether > 0
    F_all_t =Field_pad;
    shift_sources_spaceOnly
    Field_pad =F_all_t + Field_pad*Source.ReleaseTogetherFactor;
    %  clear Field
end
clear Field
if Source.Number >= 2
    Source = rmfield(Source, 'Field_padded') ; pause(0.1)
end
%% particle velocity

% vt(:,:,:) = diff(Field(:,:,:),[],3) ./ Params.delta_t   ;
% subplot(3,1,3);imagesc(squeeze(Field_static(:,:)),limits);colorbar;title(['Gstat' num2str(i)]);

%vt(:,:,:) = diff(Ut(:,:,:),[],3) ./ Params.delta_t   ;
% % % Field_pad = Field_pad(:,:,1:floor(1/Params.Fs  /Params.delta_t):size(Field_pad,3)-Params.padnumber.time);
% % % Field = Field(:,:,1:floor(1/Params.Fs  /Params.delta_t):size(Field,3));
% % % % field is cleared
if isequal(Params.Vp_Dpl,'Vp')   % decide if you want displacement or Vp as primary output
    %     vt = Field_pad(:,:,1:floor(1/Params.Fs  /Params.delta_t):size(Field_pad,3)-Params.padnumber.time);
    Ut = cumsum(Field_pad,3).*delta_t ;
    Ut = Ut(:,:,1:floor(1/Params.Fs  /Params.delta_t):size(Field_pad,3)-Params.padnumber.time);
    vt = diff(Ut,[],3) ./ (1/Params.Fs)  ;
    
elseif isequal(Params.Vp_Dpl,'Dpl')
    Ut = Field_pad(:,:,1:floor(1/Params.Fs  /Params.delta_t):size(Field_pad,3)-Params.padnumber.time);
    vt = diff(Ut,[],3) ./ (1/Params.Fs)  ;
end
'finish vt'

%% save Data
%save simulation G_fct G_fct_near G_fct_s Ut vt Source Params r Trans -v7.3
%save([Params.filename '.mat'] , 'Ut', 'Fieldterms', 'vt', 'Source', 'Params', 'r', 'Trans', '-v7.3')
if isequal(Params.saveSwitch,'on')
    save([Params.filename 'Short.mat'] ,  'Ut', 'vt', 'Source', 'Params', 'r', 'Trans', '-v7.3')
end

%% Make VIdeo

PrCTILEdpl = prctile(abs(Ut(:)),99) ;
PrCTILEvt = prctile(abs(vt(:)),99) ;
limitsdpl = [-PrCTILEdpl PrCTILEdpl];
limitsvt  = [-PrCTILEvt PrCTILEvt];
videoname =  Params.filename  ;   %Params.filename ;
subplot_pos = '1,2,'

if isequal(Params.video ,'on') == 1
    Particle_Deplacement_film(Ut,vt, Params, [] ,videoname, limitsdpl, limitsvt, subplot_pos )
end

nout = max(nargout,1) - 1;

for k = 1:nout
            if Source.ReleaseTogether > 0

    if k == 1
            varargout{k} = F_all_t;
        
    elseif k==2
        
        varargout{k} = diff(F_all_t(:,:,:),[],3) ./ (Params.delta_t)  ;
    elseif k == 3
        varargout{k} = G_fct;
    else
        display('this input does not exist')
        end
    end
end




end
