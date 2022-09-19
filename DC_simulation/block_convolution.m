function Field = block_convolution(source_ux,G_fct,block_size)
% convolutes a 1d source(time) vector with a 3d matrix
% optimizes computation time by block convoluting large matrices
% time dimension needs to be 3rd dimension
% Inputs:
% source_ux :1d vector  , G_fct - Greens fct (3d), 
% block_size - scalar value for square first 2 dimensions) 

iter_i = ceil(size(G_fct,1)/block_size) ;
iter_j = ceil(size(G_fct,2)/block_size) ;

test= reshape(source_ux,[1 1 length(source_ux)]);
Field = zeros(size(G_fct,1),size(G_fct,2),size(G_fct,3)+length(test)-1 );
for i = 1:iter_i
    
    
    if i == iter_i
        A = (i-1)*block_size +1 : size(G_fct,1);
        
    else
        A = (i-1)*block_size +1 : i*block_size;
    end
    
    for j = 1:iter_j
        if j == iter_j
            B = (j-1)*block_size +1 : size(G_fct,2);
        else
            B = (j-1)*block_size +1 : j*block_size;
        end
        tic
        Field(A,B,:) = convn(test,G_fct(A,B,:));
        toc
    end
end