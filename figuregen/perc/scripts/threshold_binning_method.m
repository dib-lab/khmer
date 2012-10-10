function [distribution,bins] = threshold_binning_method(D,T)
  n=1;
  [z]=length(D);
  b=1;
  bins(1,1)=1;
  for l=1:z
      lots=find(bins(:,1)==b);
      S=sum(D(lots));
      if(S>T)
          bins(l,1)=b;
          b=b+1;
      else
          bins(l,1)=b;
      end
  end
  distribution=get_dist_into_bins(D,bins);
end

