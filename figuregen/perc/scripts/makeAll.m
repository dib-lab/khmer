%
% This file is part of khmer, http://github.com/ged-lab/khmer/, and is
% Copyright (C) Michigan State University, 2009-2014. It is licensed under
% the three-clause BSD license; see doc/LICENSE.txt. 
% Contact: khmer-project@idyll.org
%
%make all
N=0.18:0.0001:0.1899;
for i=1:length(N)
    m=0;
    loops=9;
        for r=0:loops
            filename=strcat('dist_R_',num2str(N(i),'%10.4f'),'_12_R',num2str(r),'.txt')
            D=load(filename);
            %D=blowUpDist(D);
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
            %plot(N(i),sum((D2.^2))/sum((D1.^2))),hold on;
            %plot(N(i),sum((D2.^2))/sum((D1.^2)),'O'),hold on;
            m=m+sum((D2.^2))/sum((D1.^2));
            drawnow
        end
        %plot(N(i),m/10),hold on;
        %plot(N(i),m/10,'+red'),hold on;
        mm(i)=m/(loops+1);
 end
 P=N(find(mm==max(mm)))
 %   title(fileString);
 %saveas(gcf,'FstatisticPlot.pdf','pdf');
 %saveas(gcf,'FstatisticPlot.fig','fig');
%print -dpdf perc.pdf
exit();
