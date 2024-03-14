  
function [ar_comp,table_wp] = make_stattable(pval,wald,area2test)

table_wp =[];ar_comp =[];
for ar1 = 1 : length(area2test)-1
    for ar2 = 1 : length(area2test)
        if ar2>ar1
            n=1;pp=0;
            while pp == 0
                n = n + 1;
                pp = round(pval(ar2,ar1),n);
            end
            pp = round(pval(ar2,ar1),n+1);

            table_wp = [table_wp;{num2str(round(wald(ar2,ar1),2)) num2str(pp)}];
            ar_comp = [ar_comp ;{area2test{ar1} area2test{ar2}}];
        end
    end
end
