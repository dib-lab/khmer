%
% This file is part of khmer, http://github.com/ged-lab/khmer/, and is
% Copyright (C) Michigan State University, 2009-2014. It is licensed under
% the three-clause BSD license; see doc/LICENSE.txt. 
% Contact: khmer-project@idyll.org
%
function [N] = TBM(D,T)
[yDim xDim]=size(D);
x1=1;
I=0;
c=1;
N(1,1)=0;
N(1,2)=0;
for i=1:yDim
    n=intersect(find(D(:,1)>=x1),find(D(:,1)<=x1+I));
    S=sum(D(n,2));
    if(S>T)
        N(c,1)=(x1+(x1+I))/2;
        N(c,2)=S/(I+1);
        c=c+1;
        x1=x1+I+1;
    else
        I=I+1;
    end
end


end

