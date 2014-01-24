%
% This file is part of khmer, http://github.com/ged-lab/khmer/, and is
% Copyright (C) Michigan State University, 2009-2013. It is licensed under
% the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
%
[L4,E4]=loadPPTData('_4.txt');
[L5,E5]=loadPPTData('_5.txt');
[L6,E6]=loadPPTData('_6.txt');
[L7,E7]=loadPPTData('_7.txt');
[L8,E8]=loadPPTData('_8.txt');
[L9,E9]=loadPPTData('_9.txt');
[L10,E10]=loadPPTData('_10.txt');
[L11,E11]=loadPPTData('_11.txt');
[L12,E12]=loadPPTData('_12.txt');

p4=plot(L4(:,1),L4(:,2), 'color', 'magenta'),hold on;
p5=plot(L5(:,1),L5(:,2), 'color', 'blue'),hold on;
p6=plot(L6(:,1),L6(:,2), 'color', 'black'),hold on;
p7=plot(L7(:,1),L7(:,2), 'color', 'yellow'),hold on;
p8=plot(L8(:,1),L8(:,2), 'color', 'red'),hold on;
p9=plot(L9(:,1),L9(:,2), 'color', [1, .5, 0]),hold on;
p10=plot(L10(:,1),L10(:,2), 'color', [1, .05, .5]),hold on;
p11=plot(L11(:,1),L11(:,2), 'color', 'cyan'),hold on;
p12=plot(L12(:,1),L12(:,2), 'color', 'green'),hold on;

e4=errorbar(L4(:,1),L4(:,2),E4(:,2)),hold on;
e5=errorbar(L5(:,1),L5(:,2),E5(:,2)),hold on;
e6=errorbar(L6(:,1),L6(:,2),E6(:,2)),hold on;
e7=errorbar(L7(:,1),L7(:,2),E7(:,2)),hold on;
e8=errorbar(L8(:,1),L8(:,2),E8(:,2)),hold on;
e9=errorbar(L9(:,1),L9(:,2),E9(:,2)),hold on;
e10=errorbar(L10(:,1),L10(:,2),E10(:,2)),hold on;
e11=errorbar(L11(:,1),L11(:,2),E11(:,2)),hold on;
e12=errorbar(L12(:,1),L12(:,2),E12(:,2)),hold on;

set(e4,'color','magenta');
set(e5,'color','blue');
set(e6,'color','black');
set(e7,'color','yellow');
set(e8,'color','red');
set(e9,'color',[1, .5, 0]);
set(e10,'color',[1, .05, .5]);
set(e11,'color','cyan');
set(e12,'color','green');

set(gca,'fontsize',14);

set(get(gca,'XLabel'),'fontsize',14);
set(get(gca,'YLabel'),'fontsize',14);
set(get(gca,'XLabel'),'String','p (fraction of active vertices)');
set(get(gca,'YLabel'),'String','\theta (relative size of largest component)');
legend('K4','K5','K6','K7','K8','K9','K10','K11','K12');
%legend('K4','K5','K6','K7');
%saveas(gcf,'resultPlot.pdf','pdf');
%exit();
print -dpdf S1.pdf
exit();
