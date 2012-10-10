function [N] = blowUpDist(D)
    [yDim xDim]=size(D);
    N=sparse(0);
    for i=1:yDim
        N(D(i,1))=D(i,2);
    end
end

