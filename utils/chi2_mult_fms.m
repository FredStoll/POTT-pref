function [tbl,chi2stat,pval] = chi2_mult_fms(nN)
       
       x1= [];
       x2=[];
       chars = 'abcdefghijklmno';
       for i  = 1 : length(nN(:,1))
            x1 = [x1 ; repmat(chars(i),nN(i,2),1)];
       end
       for i  = 1 : length(nN(:,1))
            x2 = [x2 ; repmat(1,nN(i,1),1) ; repmat(2,nN(i,2)-nN(i,1),1)];
       end
       [tbl,chi2stat,pval] = crosstab(x1,x2);
       
       
       
       
       
       
       
       
       
       