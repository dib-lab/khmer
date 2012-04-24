function [L,E] = loadPPTData(trailer)
p=0.01:0.01:0.99;
for i=1:length(p)
    L(i,1)=p(i);
	E(i,1)=p(i);
    filename=strcat('R_',num2str(p(i),'%10.2f'),trailer);
    D=load(filename);
    [yd xd]=size(D);
    L(i,2:xd+1)=mean(D);
	E(i,2:xd+1)=std(D)/sqrt(yd);
end

end
