
if isequal(Rupture.rup_direction,'pos') == 1
    Source.Field_padded = padarray(Field,[0  0 padnumber_time],0,'post');
    Source.Field_padded = padarray(Field,[0  ceil(1.0*padnumber_space) 0],0,'pre');
elseif  isequal(Rupture.rup_direction,'neg') ==1
    Source.Field_padded = padarray(Field,[0  0 padnumber_time],0,'post');
    Source.Field_padded = padarray( Source.Field_padded, [0  ceil(1.0*padnumber_space) 0],0,'post');   
end
%clear Field
Field_pad =Source.Field_padded ;
% loop over all sources and shift the right amount
j=1; %iteration for the Source Point coordinate
for i = 2:rup_steps
    
    
    if ismember(i,PosOne) == 1
        %find the space shift
        [~,B] = ismember(i,PosOne) ;
        
        switch(Rupture.rup_direction)
            case{'pos'}
                Field_pad= Field_pad + circshift(Source.Field_padded, [0 B-1 i-1]);
                
            case{'neg'}
                Field_pad= Field_pad + circshift(Source.Field_padded, [0 -B+1 i-1]);
        end                               %         wave_tot_p = wave_tot_p + circshift(Source(1).wave_p, [0 B-1 i-1]);
        j = j+1;
        
%         Source.Point(j,2) = source_ind(2) +(i-1);
    end
    
    
end

