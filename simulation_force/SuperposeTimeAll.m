% run over several frequencies and
% plot
p_speed = alpha;
figure
h1 = subplot(3,2,1); legend ;
h2 = subplot(3,2,2); legend ; hold on
h3 = subplot(3,2,3); legend ;
h4 = subplot(3,2,4); legend ; hold on
h5 = subplot(3,2,5); legend ;
h6 = subplot(3,2,6); legend ; hold on


colorhandle = rand(1,3);

coloR = colorhandle(1,:);
observation_dir = obs_dir;
Source.Point(1,:) = Source_ind();
density = rho ;
[r,depth_mat,x_mat] = Grid_Fourier(depth_points,length_points,delta_x,delta_z,Source_ind) ;
Source.X0 = Delta_source_function(time(1:timesteps), r, beta,Source.Point,depth_mat,'discrete')  ;
Source.X1 = Delta_source_function(time(1:timesteps), r, p_speed,Source.Point,depth_mat,'discrete')  ;
Source.X2 = Delta_near_field(r,time,beta,p_speed) ;



wave_type = 's-wave' ;
Green_type_switch ;
wave_type = 'p-wave' ;
Green_type_switch ;
wave_type = 'near-field' ;
Green_type_switch ;
G_fct = G_fct_near + G_fct_s + G_fct_p ;
Gnear = G_fct_near; clear G_fct_near
Gs = G_fct_s ; clear G_fct_s
Gmn = G_fct ; clear G_fct ; clear G_fct_p
%shift the sources
Gmn_final = Gmn ;
Gnear_final = Gnear ;
Gs_final = Gs ;
jj = 1;
ShiftedG(jj).Gmn = Gmn ;
ShiftedG(jj).Gnear = Gnear ;
ShiftedG(jj).Gs = Gs ;
jj = jj+1 ; 
for kk = [-shift_length:-1 1:shift_length]
    ShiftedG(jj).Gmn = circshift(Gmn,[kk 0]) ;
    ShiftedG(jj).Gmn(1:shift_length,:) = nan;
    ShiftedG(jj).Gmn(end-shift_length:end,:) =  nan;
%     imagesc(abs(ShiftedG(jj).Gmn)); 
%     pause
    ShiftedG(jj).Gnear = circshift(Gnear,[kk 0]) ;
    ShiftedG(jj).Gnear(1:shift_length,:) = nan;
    ShiftedG(jj).Gnear(end-shift_length:end,:) =  nan;
    
    ShiftedG(jj).Gs = circshift(Gs,[kk 0]) ;
    ShiftedG(jj).Gs(1:shift_length,:) = nan;
    ShiftedG(jj).Gs(end-shift_length:end,:) =  nan;  
    
    Gnear_final  =  ShiftedG(jj).Gnear+Gnear_final ;    
    Gmn_final  =  ShiftedG(jj).Gmn+Gmn_final ;
    Gs_final  =  ShiftedG(jj).Gs+Gs_final ;
jj = jj+1 ;
end
GmnFreq(1).Gmn_final = Gmn_final;
GmnFreq(1).Gnear_final = Gnear_final;
GmnFreq(1).Gs_final = Gs_final;

rplot = r; rplot( Source_ind(1),Source_ind(2)) = 0;
rplotx =  rplot(Source_ind(1),:) ; rplotx(1:Source_ind(2)) = - rplotx(1:Source_ind(2))  ; 
rploty =  rplot(:,Source_ind(2)) ; rploty(1:Source_ind(1)) = - rploty(1:Source_ind(1))  ; 
% lam_s = beta/f;
% lam_p = p_speed/f;

% % figure 1
% % semilogy(rplotx , abs(Gmn_final(Source_ind(1),:)),'DisplayName',[num2str(f) ' Hz'],'Color',coloR) ;hold on ;
% semilogy(h2,rplotx , mean(abs(Gmn_final(shift_length+1:end-(shift_length+1),:)),1),...
%    'Color',coloR) ;hold on ;
% title('along x')
% % if lam_s < rplotx(end)
% %     v1 = vline(lam_s,':',[],[],[],h2) ; 
% %     v1.Color = coloR;
% % else
% %     v2 = vline(rploty(end),':',[],[],[],h2) ; 
% %     v2.Color = coloR;
% %  text(h2,  rplotx(end-5),floor(max(abs(Gmn_final(Source_ind(1),:)))),...
% %     [ '$\lambda$ = ' num2str(lam_s)])
% % end
% % 'label',num2str(f)
% 
% % figure 2
% semilogy(h4, rplotx , mean(abs(Gnear_final(shift_length+1:end-(shift_length+1),:)),1),...
%     'Color',coloR) ;hold on ;
% title(h4, 'along x')
% % if lam_s < rplotx(end) 
% %     v1 = vline(lam_s,':',[],[],[],h4) ; 
% %     v1.Color = coloR;
% % else
% %     v2 = vline(rploty(end),':',[],[],[],h4)  ; 
% %     v2.Color = coloR;
% %  text( h4, rplotx(end-5),floor(max(abs(Gnear_final(Source_ind(1),:)))),...
% %     [ '$\lambda$ = ' num2str(lam_s)])
% % end
% % 'label',num2str(f)
% 
% % figure 3
% semilogy(h6, rplotx , mean(abs(Gs_final(shift_length+1:end-(shift_length+1),:)),1),...
%    'Color',coloR) ;hold on ;
% title(h6, 'along x')
% % if lam_s < rplotx(end) 
% %     v1 = vline(lam_s,':',[],[],[],h6) ; 
% %     v1.Color = coloR;
% % else
% %     v2 = vline(rploty(end),':',[],[],[],h6)  ; 
% %     v2.Color = coloR;
% %  text( h6, rplotx(end-5),floor(max(abs(Gs(Source_ind(1),:)))),...
% %     [ '$\lambda$ = ' num2str(lam_s)])
% % end
% 
% 
% 
% 
% %figures
% imagesc(h1, depth_mat(1:end,1),x_mat(1,1:end),abs(Gmn_final)); colorbar(h1); title(h1,'Superposed Total G-fct');xlabel(h1,'x [m]');ylabel(h1,'z [m]')
% legend(h2,'show');title(h2,'Total Green fct.');xlabel(h2,'x [m]'); ylabel(h2,'Amp.')
% imagesc(h3,depth_mat(1:end,1),x_mat(1,1:end),abs(Gnear_final)); colorbar(h3); title(h3,'Superposed Near-field G-fct');xlabel(h3,'x [m]');ylabel(h3,'z [m]')
% legend(h4,'show');title(h4,'Near-field Green fct.');xlabel(h4,'x [m]'); ylabel(h4,'Amp.')
% imagesc(h5,depth_mat(1:end,1),x_mat(1,1:end),abs(Gs_final)); colorbar(h5); title(h5,'Superposed S G-fct');xlabel(h5,'x [m]');ylabel(h5,'z [m]')
% legend(h6,'show');title(h6,'S Green fct.');xlabel(h6,'x [m]'); ylabel(h6,'Amp.')
% suptitle(['Force Direction: ' obs_dir ';  Source relative to observation: ' num2str(source_dir) '$^\circ$'])
% 
% % 
% % figure
% % subplot(2,3,1)
% % imagesc(abs(Gmn))
% % subplot(2,3,2)
% % plot( rplotx , abs(Gmn(Source_ind(1),:)) ) ;hold on ;
% % title('along x')
% % if lam_s < rplotx(end)
% %     vline(lam_s)
% % else
% %      vline( rploty(end))
% %  text( rplotx(end-5),floor(max(abs(Gmn(Source_ind(1),:)))),...
% %     [ '$\lambda$ = ' num2str(lam_s)])
% % end
% % subplot(2,3,4)
% % imagesc(angle(Gmn))
% % subplot(2,3,5)
% % plot(angle(Gmn(Source_ind(1),:)))
% % subplot(2,3,3)
% % plot(rploty,abs(Gmn(:,Source_ind(2)))); hold on
% % title('along y')
% % if lam_s < rploty(end)
% %     vline(lam_s)
% % else
% %     vline( rploty(end))
% %     text( rploty(end-5),floor(max(abs(Gmn(Source_ind(2),:)))),...
% %     [ '$\lambda$ = ' num2str(lam_s)])
% % end
% % subplot(2,3,6)
% % plot(angle(Gmn(:,Source_ind(2))))
% % suptitle( ['$V_{p}$= ' num2str(p_speed) '  $V_{s}$=' num2str(beta)  ' f=' num2str(f) ' $\Delta\phi_{SourceObservation}$=' num2str(delta_angle) ])
