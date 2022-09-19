

% workspace


if isequal(Params.Vp_Dpl,'Vp')   % decide if you want displacement or Vp as primary output   
        %     vt = Field_pad(:,:,1:floor(1/Params.Fs  /Params.delta_t):size(Field_pad,3)-Params.padnumber.time);
       Field_s = block_convolution(Source.Function_deriv, G_fct_s, 60);
    Field_near = block_convolution(Source.Function_deriv, G_fct_near, 60);
    Field = block_convolution(Source.Function_deriv, G_fct, 60);
        
        Ut = cumsum(Field,3) ;
    Ut = Ut(:,:,1:floor(1/Params.Fs  /Params.delta_t):size(Field,3));
    Ut_s = cumsum(Field_s,3) ;
    Ut_s = Ut_s(:,:,1:floor(1/Params.Fs  /Params.delta_t):size(Field_s,3));
    
    Ut_near = cumsum(Field_near,3) ;
    Ut_near = Ut_near(:,:,1:floor(1/Params.Fs  /Params.delta_t):size(Field_near,3));
     
    vt = diff(Ut,[],3) ./ (1/Params.Fs)  ;
    vt_s = diff(Ut_s,[],3) ./ (1/Params.Fs)  ;
    vt_near = diff(Ut_near,[],3) ./ (1/Params.Fs)  ;    
elseif isequal(Params.Vp_Dpl,'Dpl')
    Field_s = block_convolution(Source.Function, G_fct_s, 60);
    Field_near = block_convolution(Source.Function, G_fct_near, 60);
    Field = block_convolution(Source.Function, G_fct, 60);
    Ut_s = Field_s(:,:,1:floor(1/Params.Fs  /Params.delta_t):size(Field_s,3));
    Ut_near = Field_near(:,:,1:floor(1/Params.Fs  /Params.delta_t):size(Field_near,3));
    Ut = Field(:,:,1:floor(1/Params.Fs  /Params.delta_t):size(Field,3));
    vt = diff(Ut,[],3) ./ (1/Params.Fs)  ;
        vt_s = diff(Ut_s,[],3) ./ (1/Params.Fs)  ;
           vt_near = diff(Ut_near,[],3) ./ (1/Params.Fs)  ;

end

 
figure; 
subplot(3,3,1)
imagesc(vt(:,:,round(size(vt,3)/3)))
subplot(3,3,2)
imagesc(vt_s(:,:,round(size(vt,3)/3)))
subplot(3,3,3)
imagesc(vt_near(:,:,round(size(vt,3)/3)))
subplot(3,3,4)
imagesc(Ut(:,:,round(size(Ut,3)/3)))
subplot(3,3,5)
imagesc(Ut_s(:,:,round(size(Ut,3)/3)))
subplot(3,3,6)
imagesc(Ut_near(:,:,round(size(Ut,3)/3)))
if isequal(Params.saveSwitch,'on')
    save([Params.filename 'RadPat.mat'] , '-v7.3')
end
%%
% save('RadPatForce.mat')
% Params.Vp_Dpl
% size(Field)
% workspace
% imagesc(Ut(:,:,10))
% clf
% imagesc(Ut(:,:,10))
% imagesc(Ut(:,:,20))
% imagesc(Ut(:,:,24))