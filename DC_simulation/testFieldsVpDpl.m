%%


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

Field_p = block_convolution(Source.Function_deriv, G_fct_p, 60);
Field_s = block_convolution(Source.Function_deriv, G_fct_s, 100);
Field_i_s = block_convolution(Source.Function, G_fct_I_s, 100);
Field_i_p = block_convolution(Source.Function, G_fct_I_p, 60);
Field_near_near = block_convolution(Source.Function, G_fct_near, 100);
Field_static = Source.Function(end) .* G_fct_static ;
toc

Field = Field_p + Field_s + Field_i_s(:,:,1:end-1) + Field_i_p(:,:,1:end-1) +Field_near_near(:,:,1:end-1);

%if Source.d_shift_min_points ~= 0
if length(find(Source.sourcesvector==1)) > 1
    tic
    shift_sources %gives Field_pad
    toc
    'finish shift'
    FieldPadDpl = Field_pad;

%     % Field = Field(:,:,1:floor(1/Params.Fs  /Params.delta_t):size(Field,3));
%     % % field is cleared
%     Ut = Field_pad(:,:,1:floor(1/Params.Fs  /Params.delta_t):size(Field_pad,3)-Params.padnumber.time);
%     vtinst = diff(Field_pad,3)%./delta_t;
%     vtinst = vtinst(:,:,1:floor(1/Params.Fs  /Params.delta_t):size(Field_pad,3)-Params.padnumber.time);
%     vt = diff(Ut(:,:,:),[],3) ./ (1/Params.Fs)  ;
end

%

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

Field = Field_p + Field_s + Field_i_s + Field_i_p +Field_near_near;

if length(find(Source.sourcesvector==1)) > 1
    tic
    shift_sources %gives Field_pad
    toc
    'finish shift'
    FieldPadVt = Field_pad;

%     Ut = cumtrapz(Field_pad,3).*delta_t ;
%     vtinst = (Field_pad);
%     vtinst = vtinst(:,:,1:floor(1/Params.Fs  /Params.delta_t):size(Field_pad,3)-Params.padnumber.time);
%     Ut = Ut(:,:,1:floor(1/Params.Fs  /Params.delta_t):size(Field_pad,3)-Params.padnumber.time);
%     vt = diff(Ut(:,:,:),[],3) ./ (1/Params.Fs)  ;
end
%%

figure; 
subplot(121); plot(Source.Function); hold on; plot(cumtrapz(Source.Function_deriv).*delta_t)
subplot(122); plot(Source.Function_deriv); hold on; plot(diff(Source.Function)./delta_t,'--')

%%
figure; 
subplot(121); imagesc(FieldPadDpl(:,:,500)); colorbar
subplot(122); imagesc(FieldPadVt(:,:,501)) ; colorbar
%%
figure; 
plotDpl = (FieldPadDpl) ; %./delta_t ; 
plotVp =  cumtrapz(FieldPadVt,3).*delta_t; 
subplot(121); imagesc(plotDpl(:,:,500)); colorbar
subplot(122); imagesc(plotVp(:,:,501)) ; colorbar



%%

%%


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
Function_deriv = Source.Function_deriv;
Acc = diff(Function_deriv)./delta_t;
Acc = [0 Acc] ;




Field_p = block_convolution(Source.Function_deriv, G_fct_p, 60);
Field_s = block_convolution(Source.Function_deriv, G_fct_s, 100);
Source.Function2 = [Source.Function repmat(Source.Function(end),1,300)]

Field_i_s = block_convolution(Source.Function2, G_fct_I_s, 100);
Field_i_p = block_convolution(Source.Function2, G_fct_I_p, 60);
Field_near_near = block_convolution(Source.Function2, G_fct_near, 100);

Field_p2 = block_convolution(Acc, G_fct_p, 60);
Field_s2 = block_convolution(Acc, G_fct_s, 100);
Field_i_s2 = block_convolution(Function_deriv, G_fct_I_s, 100);
Field_i_p2 = block_convolution(Function_deriv, G_fct_I_p, 60);
Field_near_near2 = block_convolution(Function_deriv, G_fct_near, 100);


%%
figure; 
subplot(121); imagesc(Field_s(:,:,500)); colorbar
subplot(122); imagesc(Field_s2(:,:,501)) ; colorbar
%%

figure; 
plotDpl = (Field_s) ; %./delta_t ; 
plotVp =  cumtrapz(Field_s2,3).*delta_t; 
subplot(121); imagesc(plotDpl(:,:,500)); colorbar
subplot(122); imagesc(plotVp(:,:,501)) ; colorbar

%%
%%
figure; 
subplot(321); imagesc(Field_s(:,:,500)); colorbar
subplot(322); imagesc(Field_s2(:,:,501)) ; colorbar
subplot(323); imagesc(Field_i_s(:,:,500)); colorbar
subplot(324); imagesc(Field_i_s2(:,:,501)) ; colorbar
subplot(325); imagesc(Field_near_near(:,:,100)); colorbar
subplot(326); imagesc(Field_near_near2(:,:,101)) ; colorbar
%%
figure; 
subplot(321); imagesc(Field_s(:,:,500)); colorbar
Field_s2Int = cumtrapz(Field_s2,3).*delta_t;
subplot(322); imagesc(Field_s2Int(:,:,501)) ; colorbar
subplot(323); imagesc(Field_i_s(:,:,500),[-10 10]); colorbar
Field_i_s2Int = cumtrapz(Field_i_s2,3).*delta_t;
subplot(324); imagesc(Field_i_s2Int(:,:,501),[-10,10]) ; colorbar
subplot(325); imagesc(Field_near_near(:,:,100),[-10 10]); colorbar
Field_near_near2Int = cumtrapz(Field_near_near2,3).*delta_t;
subplot(326); imagesc(Field_near_near2Int(:,:,101),[-10 10]) ; colorbar

%%

figure; 
plotDpl = (Field_s) ; %./delta_t ; 
plotVp =  cumtrapz(Field_s2,3).*delta_t; 
subplot(121); imagesc(plotDpl(:,:,500)); colorbar
subplot(122); imagesc(plotVp(:,:,501)) ; colorbar

%%
%%
figure; 
subplot(121); imagesc(Field_s(:,:,500)); colorbar
subplot(122); imagesc(Field_s2(:,:,501)) ; colorbar
%%

figure; 
plotDpl = (Field_s) ; %./delta_t ; 
plotVp =  cumtrapz(Field_s2,3).*delta_t; 
subplot(121); imagesc(plotDpl(:,:,500)); colorbar
subplot(122); imagesc(plotVp(:,:,501)) ; colorbar
%%
figure; 
plotDpl = (Field_near_near) ; %./delta_t ; 
plotVp =  cumtrapz(Field_near_near2,3).*delta_t; 
subplot(121); imagesc(plotDpl(:,:,500)); colorbar
subplot(122); imagesc(plotVp(:,:,501)) ; colorbar

