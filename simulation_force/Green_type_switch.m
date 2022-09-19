switch(wave_type)
    case{'s-wave'}
        G_fct_s =  Green_function(observation_dir,Params.observation_plane,...
            beta,Source.X0,timesteps,...
            r,Source.Point(1,:),depth_mat,density,wave_type,time,Source.Direction);
        %%convolve
    case{ 'p-wave' }
        %P-field
        G_fct_p = Green_function(observation_dir,Params.observation_plane,...
            p_speed,Source.X1,timesteps,...
            r,Source.Point(1,:),depth_mat,density,wave_type,time,Source.Direction);
        %         Source(1).wave_p = padarray(wave_p,[0  padnumber_space padnumber_time],0,'post');
        %         wave_tot_p = Source(1).wave_p ; toc
    case{ 'near-field' }
  %      Source.X2 = zeros(size(r,1), size(r,2), length(time));
%%%        for t = 1:length(time)-1
  %          X2 = zeros(size(r,1), size(r,2)) ;
  %          testplus_s = zeros(size(r,1), size(r,2));
  %          testplus_p = zeros(size(r,1), size(r,2));
  %
  %          testplus_s(:,:) = time(t).*ones(size(r)) - r(:,:)./beta;
  %          testplus_p(:,:) = time(t).*ones(size(r)) - r(:,:)./p_speed;
  %          Signp = sign(testplus_p) ;
  %          [B] = find(Signp == 1);
  %          Signs = sign(testplus_s) ;
  %          [A] = find(Signs == -1) ;
  %%          [C] = intersect(A,B) ;
    %        X2(C) = 1 ;
    %        Source.X2(:,:,t) = X2 ;
    %    end
      %  Source.Direction
      %  G_fct_near = Green_function(observation_dir,beta,Source.X2,timesteps,r,Source.Point(1,:),depth_mat,density,wave_type,time);
        G_fct_near = Green_function(observation_dir,Params.observation_plane,...
            beta,Source.X2,timesteps,...
            r,Source.Point(1,:),depth_mat,density,wave_type,time,Source.Direction);
end
