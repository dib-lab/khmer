%
% This file is part of khmer, http://github.com/ged-lab/khmer/, and is
% Copyright (C) Michigan State University, 2009-2013. It is licensed under
% the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
%
function [P,N,mm] = findMaxF(fileString)
    N=0.10:0.01:0.59;

    for i=1:length(N)
        m=0;
        for r=0:9
            D=load(strcat('dist_R_',num2str(N(i),'%10.2f'),fileString,'R',num2str(r),'.txt'));
            x=D(:,1);
            y=D(:,2);
            [p2,S2,mu2] = polyfit(log(x),log(y),2);
            [p1,S1,mu1] = polyfit(log(x),log(y),1);
            pv2=polyval(p2,log(x));
            pv1=polyval(p1,log(x));
            D1=pv2-log(y);
            D2=pv1-log(y);
            plot(N(i),sum(abs(D2))/sum(abs(D1)),'O'),hold on;
            m=m+sum(abs(D2))/sum(abs(D1));
            drawnow
        end
        plot(N(i),m/10,'+red');
        mm(i)=m/10;
    end
    P=N(find(mm==max(mm)));
end
