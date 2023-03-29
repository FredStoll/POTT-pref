function [lme,model_final] = model_comparison(modeldata_su,models_form,binom)

if binom
    lme_sess = fitglme(modeldata_su,models_form{1},'Distribution','Binomial');
    lme_sess2 = fitglme(modeldata_su,models_form{2},'Distribution','Binomial');
else
    lme_sess = fitglme(modeldata_su,models_form{1});
    lme_sess2 = fitglme(modeldata_su,models_form{2});    
end
%- pick the lowest model
nope = false;
try res = compare(lme_sess2,lme_sess);
catch nope = true;
end
if nope
    lme = lme_sess; model_final = models_form{1};
else
    if res.pValue<0.05
        lme = lme_sess; model_final = models_form{1};
    else
        lme = lme_sess2; model_final = models_form{2};
    end
end
