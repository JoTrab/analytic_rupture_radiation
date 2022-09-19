function  Field_near_near = continuous_displacement(Field_near_near)
%continuous_displacement takes the maximum or minimum displacement from the
%convolution and continues it in time until the end of the simulation
%   

[Max_near,Ind]= max(abs(Field_near_near),[],3);
for i =1:size(Field_near_near,1)
    for j = 1:size(Field_near_near,2)
        Sign_max_near(i,j) = sign(Field_near_near(i,j,Ind(i,j)));
        Field_near_near(i,j,Ind(i,j):end) = repmat(Sign_max_near(i,j) .* Max_near(i,j),[1,1,size(Field_near_near,3) - (Ind(i,j)-1)]);
       
    end
end

end

