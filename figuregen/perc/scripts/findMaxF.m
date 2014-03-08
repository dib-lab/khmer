%
% This file is part of khmer, http://github.com/ged-lab/khmer/, and is
% Copyright (C) Michigan State University, 2009-2013. It is licensed under
% the three-clause BSD license; see doc/LICENSE.txt. 
% Contact: khmer-project@idyll.org
%
function [P,N,mm] = findMaxF(fileString,doMax)
    N=0.010:0.001:doMax;
    %N=0.001:0.001:0.300;
    clf('reset');
    for i=1:length(N)
        m=0;
        loops=9;
        for r=0:loops
            D=load(strcat('data/dist_R_',num2str(N(i),'%10.3f'),fileString,'R',num2str(r),'.txt'));
            if(length(D)<3)
                D1=1;
                D2=0;
            else
                [U]=TBM(D,100);
                x=U(:,1);
                y=U(:,2);
                [p2,S2,mu2] = polyfit(log(x),log(y),2);
                [p1,S1,mu1] = polyfit(log(x),log(y),1);
                pv2=polyval(p2,log(x));
                pv1=polyval(p1,log(x));
                D1=pv2-log(y);
                D2=pv1-log(y);
            end
            subplot(2,1,1),plot(N(i),sum((D2.^2))/sum((D1.^2)),'O'),hold on;
            m=m+sum((D2.^2))/sum((D1.^2));
            drawnow
        end
        subplot(2,1,2),plot(N(i),m/10,'+red'),hold on;
        mm(i)=m/(loops+1);
    end
    P=N(find(mm==max(mm)));
    title(fileString);
    saveas(gcf,strcat('figure_',fileString,'.pdf'),'pdf');
    saveas(gcf,strcat('figure_',fileString,'.fig'),'fig');
end
