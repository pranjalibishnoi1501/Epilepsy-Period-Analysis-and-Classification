%%%%% Learn graph from the data (distance method) %%%%%%

function [A,W,Top] = learnGraph(data,sigma,N)

    W = zeros(N,N);

    %% Set weights
    for i = 1:N
        for j = 1:N
            if(i~=j)
                %W(i,j) = exp(-(norm(data(:,:,i)-data(:,:,j),2))/(2*(sigma^2))); % Vector distances
                W(i,j) = corr2(data(:,:,i), data(:,:,j));
            end
        end
    end
    
    %% Threshold
    A = W;
    [A,W,Top] = thres(A,N,N);
    W = normAdj(W);

end