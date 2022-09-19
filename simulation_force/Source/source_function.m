function X0 = source_function(time, r, beta, source)
% SOURCE function , gives the source effect over time for a 2d plane -->
% does the discretization. says which signal has arrived or has not
% arrived. Meaning that if time - r(:,:)/speed is close to 0 it is valid (?)
% USAGE:
%       X0 = source_function(time, r, beta, source)
%
% INPUTS:
%       time   - input time array [s]
%       r  - distance vector field from source (source - receiver vector) [m]
%       beta  - wave speed  [m/s]
%       source  - source characteristics. needed at least:source
%


X0 = zeros(size(r,1), size(r,2), length(time));





parfor t = 1:length(time)-1
    %     if t <= length(source)-1



    X_iter =   zeros(size(r,1), size(r,2));
    testminus = zeros(size(r,1), size(r,2));
    testplus = zeros(size(r,1), size(r,2));


    %to test if the time of flight at this timestep is negative,0
    %or positive
    %
    testplus = time(t).*ones(size(r)) - r(:,:)./beta;
    Sign1 = sign(testplus) ;
    [B] = find(Sign1 == 1);

    if t == 1

        X_iter(B) = 1 ;

    else
        testminus =  time(t-1).*ones(size(r)) - r(:,:)./beta;


        Sign_1 = sign(testminus) ;
        [A] = find(Sign_1 == -1);
        [IND] =  intersect(A,B);  %find the points where the sign ...
       %changes and thus timepoint equals wave traveltime


       X_iter(IND) = 1 %one because of delta    %source(s_iter); %at these points the sourcesignal
        %of s_iter arrives for the timestep t in time


        %         imagesc(squeeze(X0(:,:,t-1)))
        %         pause
        %                 if s_iter > 100
        %                                 imagesc(X_iter,[-1 1])
        %                                   pause(0.0001)
        %                 end
    end
                [C] = find(Sign1 == 0);

        if isempty(C) == 0
            X_iter(C) = 1 %one because of delta %source(s_iter);
        end

    X0(:,:,t) = X_iter;


    %     else this part for separate calculation can be found in backupfile
    %
    %take the mean or sum of all the source times that influence at this point of time
    %        imagesc(X0(:,:,t),[-0.001 0.001]); colorbar; title(t);pause(0.000001);
    %     pause
end
