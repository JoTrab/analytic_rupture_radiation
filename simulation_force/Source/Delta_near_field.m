function X0 = Delta_near_field(r,time,beta,p_speed)

X0 = single(zeros(size(r,1), size(r,2),1));
X0 = repmat(X0,1,1,length(time));
% X0 = zeros(size(r,1), size(r,2), length(time));
for t = 1:length(time)-1
    X2 = zeros(size(r,1), size(r,2)) ;
    testplus_s = zeros(size(r,1), size(r,2));
    testplus_p = zeros(size(r,1), size(r,2));
    
    testplus_s(:,:) = time(t).*ones(size(r)) - r(:,:)./beta;
    testplus_p(:,:) = time(t).*ones(size(r)) - r(:,:)./p_speed;
    Signp = sign(testplus_p) ;
    [B] = find(Signp == 1);
    Signs = sign(testplus_s) ;
    [A] = find(Signs == -1) ;
    
    [C] = intersect(A,B) ;
    X2(C) = 1 ;
    X0(:,:,t) = X2 ;
end
end