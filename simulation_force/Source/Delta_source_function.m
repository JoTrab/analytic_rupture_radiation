function X0 = Delta_source_function(time, r, beta,Source_Point,depth_mat,calc_type)
% SOURCE function , gives the source effect over time for a 2d plane -->
% does the discretization. says which signal has arrived or has not
% arrived. If the sign changes from one timestep to next
%time - r(:,:)/speed is close to 0 it is valid --> delta is 1
% USAGE:
%       X0 = source_function(time, r, beta, source)
%
% INPUTS:
%       time   - input time array [s]
%       r  - distance vector field from source (source - receiver vector) [m]
%       beta  - wave speed  [m/s]
%       source  - source characteristics. needed at least:source
%
%%
dt = abs(time(2)-time(1)) ;
X0 = single(zeros(size(r,1), size(r,2),1));
X0 = repmat(X0,1,1,length(time));
% X0 = zeros(size(r,1), size(r,2), length(time));
switch(calc_type)
    case{'discrete'}      
        for t = 1:length(time)-1
            %     if t <= length(source)-1
            %%
            X_iter =   zeros(size(r,1), size(r,2));
            testminus = zeros(size(r,1), size(r,2));
            test = zeros(size(r,1), size(r,2));

            %to test if the time of flight at this timestep is negative,0
            %or positive
            
            test = time(t).*ones(size(r)) - r./beta;
            Sign = sign(test) ; B = find(Sign == 1);
            
            
            if t == 1
                X_iter(B) = 1./dt ;
                
            else
                testminus = time(t-1).*ones(size(r)) - r./beta;
                SignMinus = sign(testminus) ; [Bminus] = find(SignMinus == -1);
                
                [IND] =  intersect(B,Bminus);  %find the points where the sign ...
                
                %if no changes with the neighbours in either of the tests --> no add_d --> add_d = 0
                
                %this needs to be scaled by the
                %distance covered in the test
                
                X_iter(IND) = 1./dt ; %one because of delta    %source(s_iter); %at these points the sourcesignal
                %of s_iter arrives for the timestep t in time
                
                %%
            end
            %     [C] = find(SignIn == 0);
            %
            %     if isempty(C) == 0
            %         X_iter(C) = delta_factor/dt ;%one because of delta %source(s_iter);
            %         %delta function is scaled by dt
            %         %and dx - space. if space infinitely
            %         %small: 1/dt; if time inf small: 1
            %     end
            
            X0(:,:,t) = X_iter;
   
        end
   
    case{'smooth'}
   
        x_distance = r(Source_Point(1),:) ; x_distance = repmat(x_distance,size(r,1),1) ;
        
        dx = abs(x_distance(1,2)-x_distance(1,1));
        dz = abs(depth_mat(1,1)-depth_mat(2,1)) ;
        
        
        % delta_factor should be around 0.2 to work properly    ==  static;
        C = (1:size(r,1)*size(r,2)) ;
        [ic,icd] = ixneighbors(r,C); % (gives linear indexes of those values (ic); and their neighbors  (icd) for all values
        ic_icd = [ic,icd]; %create a matrix containing indices and indices of neighbours
        ic_icd = sortrows(ic_icd,1); %sort according to indices and get the necessary neighbors
        ic= ic_icd(:,1);
        icd= ic_icd(:,2);
        % ic and icd are now sorted
        ic_dis_neg = find(x_distance(ic) ~= x_distance(icd) ); %gives indices of where dx is ...
        ic_dis_pos = find(x_distance(ic) == x_distance(icd) );
        %unequal for the value and its neighbors
        ic_depth_neg = find(depth_mat(ic) ~= depth_mat(icd) ); %should work
        ic_depth_pos = find(depth_mat(ic) == depth_mat(icd) ); %should work
        
        diag = intersect(ic_depth_neg,ic_dis_neg);
        ic_diag = ic(diag); %returns  values which in the two vectors
        icd_diag = icd(diag) ;
        %which are the indices of ic where both x and z are different
        %---> returns diagonal members
        z = intersect(ic_depth_neg,ic_dis_pos) ;
        ic_z = ic(z); %returns where x is equal and z is different
        icd_z = icd(z) ;
        x = intersect(ic_depth_pos,ic_dis_neg);
        ic_x = ic(x); %returns where z is equal and x is different
        icd_x = icd(x);
        
        r_diag = sqrt(2*(r(ic_diag).^2 + r(icd_diag).^2) ./ 4) ;
        r_z = sqrt(2*(r(ic_z).^2 + r(icd_z).^2) ./ 4) ;    %length of seitenhalbierende in z --> new r
        r_x = sqrt(2*(r(ic_x).^2 + r(icd_x).^2) ./ 4)   ;  %length of seitenhalbierende in x --> new r
        
        diag_half_aligned = 0.5*sqrt((dx^2 + dz^2));
        %assign those that are not triangles but on one line (upper and lower corner)
        halfdistances = 0.5*((r(ic_diag) - r(icd_diag))) ;
        ind = find(abs(abs(halfdistances)-diag_half_aligned) <= 1*10^-15); %those are the corner values
        %     r_diag(ind) = r(icd_diag(ind)) - diag_half_aligned; %testing the snaller values is sufficient
        
        neg = find(halfdistances(ind) < 0)  ;
        pos = find(halfdistances(ind) > 0) ;
        r_diag(ind(pos)) = r(ic_diag(ind(pos))) - diag_half_aligned; % point is further than neighbor
        r_diag(ind(neg)) = r(ic_diag(ind(neg))) + diag_half_aligned; % point is closer than neighbor
        
        ic_r = [[ic_x;ic_z;ic_diag,],[r_x;r_z;r_diag]];
        ic_r = sortrows(ic_r,1);
        
        MaxInd = max(ic);
        [ic_unique,icu_da icu_dx] = unique(ic);
        %the following to lines do the same vectorized as the loop below in
        %comments
        rmax = accumarray(icu_dx,ic_r(:,2),[],@max);
        rmin = accumarray(icu_dx,ic_r(:,2),[],@min);
        tof= (rmax-rmin)./beta ;
        delta_factor = dt./(tof);  %delta function is scaled by dt (timestep) since delta is inf
        
        %     for jj = 1:length(ic_unique)
        %         ii = ic_unique(jj);
        %         ind = find(ic_r(:,1) == ii) ;
        %         rmax(jj) = max(ic_r(ind,2)) ;
        %         rmin(jj) = min(ic_r(ind,2)) ;
        %     end%will have size of ic
        %%
        for t = 1:length(time)-1
            %     if t <= length(source)-1
            
            
            %%
            X_iter =   zeros(size(r,1), size(r,2));
            testminus = zeros(size(r,1), size(r,2));
            test_neigh_min = zeros(size(r,1), size(r,2));
            test_neigh_max = zeros(size(r,1), size(r,2));
            
            
            
            %to test if the time of flight at this timestep is negative,0
            %or positive
            
            test_neigh_min = time(t).*ones(size(rmin)) - rmin./beta;
            SignNeigh_min = sign(test_neigh_min) ; [BNeigh_min] = find(SignNeigh_min == 1);
            
            test_neigh_max = time(t).*ones(size(rmax)) - rmax./beta;
            SignNeigh_max = sign(test_neigh_max) ; [BNeigh_max] = find(SignNeigh_max == -1);
            
            
            %     t_o_f_cell = abs(max(add_d)) / beta ; %time of flight over 1 gridcell (diagonal)
            %     BNeigh = union(ic_x(BNeigh_x),union(ic_z(BNeigh_z),ic_diag(BNeigh_diag)));
            Btot = (intersect(ic_unique(BNeigh_min),ic_unique(BNeigh_max))) ;
            %%
            %          t_o_f_cell = diag_half_aligned / beta ; %time of flight of arival over 1 gridcell (diagonal)
            
            
            if t == 1
                X_iter(Btot) = delta_factor(Btot)./dt ;
                
            else
                test_neigh_max = time(t-1).*ones(size(rmax)) - rmax./beta;
                SignNeigh_max = sign(test_neigh_max) ; [BNeigh_max] = find(SignNeigh_max == -1);
                
                [IND] =  intersect(Btot,ic_unique(BNeigh_max));  %find the points where the sign ...
                
                %if no changes with the neighbours in either of the tests --> no add_d --> add_d = 0
                
                %this needs to be scaled by the
                %distance covered in the test
                
                X_iter(IND) = delta_factor(IND)./dt ; %one because of delta    %source(s_iter); %at these points the sourcesignal
                %of s_iter arrives for the timestep t in time
                
                %%
            end
            %     [C] = find(SignIn == 0);
            %
            %     if isempty(C) == 0
            %         X_iter(C) = delta_factor/dt ;%one because of delta %source(s_iter);
            %         %delta function is scaled by dt
            %         %and dx - space. if space infinitely
            %         %small: 1/dt; if time inf small: 1
            %     end
            
            X0(:,:,t) = X_iter;   
        end
end
