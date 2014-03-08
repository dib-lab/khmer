%
% This file is part of khmer, http://github.com/ged-lab/khmer/, and is
% Copyright (C) Michigan State University, 2009-2013. It is licensed under
% the three-clause BSD license; see doc/LICENSE.txt. 
% Contact: khmer-project@idyll.org
%
function [Ps] = findAllP()
fileExt={'_P1.0_K12_','_P0.8_K12_','_P0.6_K12_','_P0.4_K12_','_P0.2_K12_','_P0.1_K12_','_P0.05_K12_','_P0.01_K12_'};
X=[1.0,0.8,0.6,0.4,0.2,0.1,0.05,0.01];
maxes=[0.3,0.2,0.2,0.1,0.1,0.05,0.01];
%for j=1:length(fileExt)
for j=3:3
    [p,N,mm]=findMaxF(char(fileExt(j)),maxes(j));
    %subplot(3,1,3),plot(X(j),p,'+'),hold on;
    save(strcat('N',char(fileExt(j)),'.mat'),'N');
    save(strcat('mm',char(fileExt(j)),'.mat'),'mm');
    save(strcat('p',char(fileExt(j)),'.mat'),'p');
    Ps(j,1)=X(j);
    Ps(j,2)=p;
end

end

