%
% This file is part of khmer, http://github.com/ged-lab/khmer/, and is
% Copyright (C) Michigan State University, 2009-2013. It is licensed under
% the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
%
function [N] = blowUpDist(D)
    [yDim xDim]=size(D);
    N=sparse(0);
    for i=1:yDim
        N(D(i,1))=D(i,2);
    end
end

