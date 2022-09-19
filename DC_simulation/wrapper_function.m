
function [Ut, vt, Fieldterms ,  Source, Params, Rupture, varargout] =...
    wrapper_function(Trans, Source, Params,  Rupture)
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

nout = max(nargout,1) - 1;

rup_speed_start_factor = Rupture.rup_speed_start_factor ;
rup_speed_end_factor = Rupture.rup_speed_end_factor ; %sqrt(2)  % sqrt(2) ;  %sqrt(2)
% max_rup_steps = Rupture.max_rup_steps ; %maximum time_points wished for calculation
% add_rupEnd_source_nr = Rupture.add_rupEnd_source_nr  ; %add number of sources / points/ timesteps that work at supershear speed
rup_length_points = round(Rupture.length_total / ((Trans.spacingMm/1000) ));
Rupture.rup_length_points = rup_length_points;
Rupture.add_rupEnd_source_nr = round(Rupture.final_speed  / ((Trans.spacingMm/1000) )); %add number of sources / points/ timesteps in meter that work at final speed
Rupture.add_rupStart_source_nr = round(Rupture.start_speed  / ((Trans.spacingMm/1000) )); %add number of sources / points/ timesteps in meter that work at final speed
Rupture.accelerating_length_points = rup_length_points - ...
    Rupture.add_rupStart_source_nr   - Rupture.add_rupEnd_source_nr;

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
    Field_pad = []  ;
    Fieldterms = []  ;
    vt = []  ;
    F_all_t = [] ;
    Ut = [] ;
    for k = 1:nout
        if k == 1
            varargout{k} = F_all_t;
        elseif k==2
            
            varargout{k} = diff(F_all_t(:,:,:),[],3) ./ (Params.delta_t)  ;
        else
            display('this input does not exist')
        end
    end
    return
end
%depth mat is y in the case of x-y
Source.X0 = Delta_source_function(time(1:timesteps), r, beta,Source.Point,depth_mat,calc_type)  ;
Source.X1 = Delta_source_function(time(1:timesteps), r, p_speed,Source.Point,depth_mat,calc_type)  ;
Source.X2 = Delta_near_field(r,time,beta,p_speed) ;

wave_type = 's-wave' ;
Green_type_switch ;
wave_type = 'p-wave' ;
Green_type_switch ;
wave_type = 'near-field' ;
Green_type_switch ;
wave_type = 'inter_p' ;
Green_type_switch ;
wave_type = 'inter_s' ;
Green_type_switch ;
wave_type = 'static' ;
Green_type_switch ;
G_fct_far = G_fct_s + G_fct_p ;
G_fct_i_near = G_fct_near + G_fct_I_s  + G_fct_I_p ;
Fieldterms.G_all = G_fct_far+ G_fct_i_near ;
Fieldterms.G_fct_far = G_fct_far;
Fieldterms.G_fct_i_near = G_fct_i_near
Fieldterms.G_fct_near = G_fct_near ;
Fieldterms.G_fct_I_s = G_fct_I_s ;
Fieldterms.G_fct_I_p = G_fct_I_p ;
Fieldterms.G_fct_p = G_fct_p ;
Fieldterms.G_fct_s = G_fct_s ;
% Source = rmfield(Source,{'X0','X1','X2'});
if exist('Params.returnSwitch') ==1
    if isequal(Params.returnSwitch ,'Green') == 1
        Field_pad = []  ;
        vt = []  ;
        F_all_t = [] ;
        Ut = [] ;
        for k = 1:nout
            if k == 1
                varargout{k} = F_all_t;
            elseif k==2
                
                varargout{k} = diff(F_all_t(:,:,:),[],3) ./ (Params.delta_t)  ;
            else
                display('this input does not exist')
            end
        end
        
        return
    end
end

Source = rmfield(Source, {'X0','X1', 'X2'}) ;


if isequal(Params.Vp_Dpl,'Dpl')
    tic
    Field_p = block_convolution(Source.Function_deriv, G_fct_p, 60);
    Field_s = block_convolution(Source.Function_deriv, G_fct_s, 100);
    Field_i_s = block_convolution(Source.Function, G_fct_I_s, 100);
    Field_i_p = block_convolution(Source.Function, G_fct_I_p, 60);
    Field_near_near = block_convolution(Source.Function, G_fct_near, 100);
    Field_static = Source.Function(end) .* G_fct_static ;
    toc
elseif isequal(Params.Vp_Dpl,'Vp')
    
    Function_deriv = Source.Function_deriv;
    % Function_deriv = [0 Function_deriv] NOT in case of TRIANGLE
    if isequal(Source.type_switch,'triangle') == 0
        Function_deriv = [0 Function_deriv] ;
        %triangle already has the right length?
    end
    Acc = diff(Function_deriv)./delta_t;
    Acc = [0 Acc] ;
    tic
    Field_p = block_convolution(Acc, G_fct_p, 60);
    Field_s = block_convolution(Acc, G_fct_s, 100);
    Field_i_s = block_convolution(Function_deriv, G_fct_I_s, 100);
    Field_i_p = block_convolution(Function_deriv, G_fct_I_p, 60);
    Field_near_near = block_convolution(Function_deriv, G_fct_near, 100);
    Field_static = Function_deriv(end) .* G_fct_static ;
    toc
end

clear G_fct_far G_fct_i_nearG_fct_near G_fct_I_s G_fct_I_p G_fct_p G_fct_s 


%Vp
if isequal(Params.Vp_Dpl,'Vp')
    Field = Field_p + Field_s + Field_i_s + Field_i_p +Field_near_near;
%Dpl
elseif isequal(Params.Vp_Dpl,'Dpl')
    Field = Field_p + Field_s + Field_i_s(:,:,1:end-1) + Field_i_p(:,:,1:end-1) +Field_near_near(:,:,1:end-1);
end
clear Field_p  Field_s  Field_i_s Field_i_p Field_near_near



%%


if exist('Params.returnSwitch') ==1
    if isequal(Params.returnSwitch ,'Field') == 1
        Ut = [];
        vt = [] ;
        F_all_t = [] ;
        for k = 1:nout
            if k == 1
                varargout{k} = F_all_t;
            elseif k==2
                
                varargout{k} = diff(F_all_t(:,:,:),[],3) ./ (Params.delta_t)  ;
            else
                display('this input does not exist')
            end
        end
        return
    end
end

% Field_2 = Field(:,:,1:floor(1/Params.Fs  /Params.delta_t):size(Field,3)); %low temporal reolution
% vt_Field2 = diff(Field_2(:,:,:),[],3)./ (1/Params.Fs)  ;
% vt_Field = diff(Field_high(:,:,:),[],3)./delta_t ;%./ (1/Params.Fs)  ;
'finish Field'

% delete green functions and fieldterms

%%
%if Source.d_shift_min_points ~= 0
if length(find(Source.sourcesvector==1)) > 1
    tic
    shift_sources %gives Field_pad
    toc
    'finish shift'
    F_all_t =Field_pad;
    % Field = Field(:,:,1:floor(1/Params.Fs  /Params.delta_t):size(Field,3));
    % % field is cleared
    if isequal(Params.Vp_Dpl,'Dpl')
        Ut = Field_pad(:,:,1:floor(1/Params.Fs  /Params.delta_t):size(Field_pad,3)-Params.padnumber.time);
%         vtinst = diff(Field_pad,3)%./delta_t;
%         vtinst = vtinst(:,:,1:floor(1/Params.Fs  /Params.delta_t):size(Field_pad,3)-Params.padnumber.time);
        vt = diff(Ut(:,:,:),[],3) ./ (1/Params.Fs)  ;
    elseif isequal(Params.Vp_Dpl,'Vp')
        Ut = cumtrapz(Field_pad,3).*delta_t ;
%         vtinst = (Field_pad);
%         vtinst = vtinst(:,:,1:floor(1/Params.Fs  /Params.delta_t):size(Field_pad,3)-Params.padnumber.time);   
        Ut = Ut(:,:,1:floor(1/Params.Fs  /Params.delta_t):size(Field_pad,3)-Params.padnumber.time);
        vt = diff(Ut(:,:,:),[],3) ./ (1/Params.Fs)  ;
    end
    
else
    F_all_t = Field ;
    display('One Source point only')
    if isequal(Params.Vp_Dpl,'Dpl')
        Field = Field(:,:,1:floor(1/Params.Fs  /Params.delta_t):size(Field,3));
        vt = diff(Field(:,:,:),[],3) ./ (1/Params.Fs)  ;
        Ut = Field ;
    elseif isequal(Params.Vp_Dpl,'Vp')
        Ut = cumtrapz(Field,3).*delta_t ;  ;
        Ut = Ut(:,:,1:floor(1/Params.Fs  /Params.delta_t):size(Field,3));
        vt = diff(Ut(:,:,:),[],3) ./ (1/Params.Fs)  ;
    end
end

'finish vt'



%% particle velocity

% vt(:,:,:) = diff(Field(:,:,:),[],3) ./ Params.delta_t   ;
% subplot(3,1,3);imagesc(squeeze(Field_static(:,:)),limits);colorbar;title(['Gstat' num2str(i)]);

%% save Data
%rmfield(Source,'Field_padded') ;
%save simulation G_fct G_fct_near G_fct_s Field_pad vt Source Params r Trans -v7.3
%save([Params.filename '.mat'] , 'Field_pad', 'Fieldterms', 'vt', 'Source', 'Params', 'r', 'Trans', '-v7.3')
if exist('Field_padded')==1
Source = rmfield(Source, 'Field_padded') ; pause(0.1)
end
save([Params.filename 'Short.mat'] ,  'vt','Ut','Source','Params', 'r', 'Trans', '-v7.3')

%% Make VIdeo
PrCTILEdpl = prctile(abs(Ut(:)),99) ;
PrCTILEvt = prctile(abs(vt(:)),99) ;


limitsdpl = [-PrCTILEdpl PrCTILEdpl];
limitsvt  = [-PrCTILEvt PrCTILEvt];
videoname =  Params.filename ;   %Params.filename ;
subplot_pos = '1,2,'


if isequal(Params.video ,'on') == 1
    Particle_Deplacement_film(Ut,vt, Params, [] ,videoname, limitsdpl, limitsvt, subplot_pos )
end


for k = 1:nout
    if k == 1
        varargout{k} = F_all_t;
    elseif k==2
        
        varargout{k} = diff(F_all_t(:,:,:),[],3) ./ (Params.delta_t)  ;
    else
        display('this input does not exist')
    end
end




end
