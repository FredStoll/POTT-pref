%- ft/rt same are wrong, as they are not normalized before cutting the
%session!!!!!

function [ALL,param] = behav_model_bins_deval(filename,perm,showfig)

%- load behav variable
TT = load(filename);

deval_bin=zeros(1,length(TT.t_FPfix));
deval_bin(1,1:find(diff(TT.t_FPfix)>500))=1;
deval_bin(1,find(diff(TT.t_FPfix)>500)+1:end)=2;


    when = @(arg,cd) TT.TrialType(strcmp(TT.TrialType_header,arg),:)==cd ;
    diff_juice = (when('I_juiceL',1) & when('I_juiceR',2)) | (when('I_juiceL',2) & when('I_juiceR',1)); %- take only the different juice trials
    restr = when('task',2) & when('brks',0) & diff_juice; %- take only INS trials (task == 2), completed (brks == 0) and diff juice

rt = TT.t_respfix-TT.t_resp(1,:);
ft = TT.t_FPfix_first-TT.t_FPon(1,:);
RT_minmax = [min(log(rt(restr))) max(log(rt(restr))) ];
FT_minmax = [min(log(ft(restr))) max(log(ft(restr))) ];

norm_RT = @(data) -1+((data-RT_minmax(1))*2)/(RT_minmax(2)-RT_minmax(1)) ;
norm_FT = @(data) -1+((data-FT_minmax(1))*2)/(FT_minmax(2)-FT_minmax(1)) ;


for dd = 1 : 2
    clearvars -except deval_bin dd filename perm showfig ALL param TT norm_RT norm_FT
    TrialType_header = TT.TrialType_header;
    
        param2keep = {'TrialType' 't_FPon' 't_FPfix' 't_FPfix_first' 't_respfix' 't_resp' 't_brks' 't_fboff' 't_fbon' 't_resp' 't_rew' 't_stim' 'task'};
        for pp = 1 : length(param2keep)
            eval([param2keep{pp} '= TT.' param2keep{pp} '(:,deval_bin==dd);']);
        end


    %- some utils
    take = @(arg) strcmp(TrialType_header,arg);
    norm = @(data) -1+((data-min(data))*2)/(max(data)-min(data)) ;
    when = @(arg,cd) TrialType(strcmp(TrialType_header,arg),:)==cd ;

    %- permutations
    if perm

        pav = when('task',1); %- leave the PAV trial alone
        ntr4perm = sum(~pav);
        tr_order = randperm(ntr4perm); %- perm the others (ins + dyn)

        param2perm = {'TrialType' 't_FPon' 't_FPfix_first' 't_respfix' 't_resp'};
        for pp = 1 : length(param2perm)
            eval(['dumm = ' param2perm{pp} '(:,~pav);']);
            eval([param2perm{pp} '=[' param2perm{pp} '(:,pav),dumm(:,tr_order)];']);
        end

        %- reboot the function to send the new TrialType matrix!
        when = @(arg,cd) TrialType(strcmp(TrialType_header,arg),:)==cd ;
    end


    %- trials to consider
    diff_juice = (when('I_juiceL',1) & when('I_juiceR',2)) | (when('I_juiceL',2) & when('I_juiceR',1)); %- take only the different juice trials
    restr = when('task',2) & when('brks',0) & diff_juice; %- take only INS trials (task == 2), completed (brks == 0) and diff juice

    %- variable to predict: ChooseJ1 - 1 for J1 / 0 for J2
    choice = abs(TrialType(take('I_chosenjuice'),restr)-2);

    if length(choice) > 100  %- trash sessions with less than 100 trials diff juice

        %% Step 1 - GLM on Choices in different juice trials
        %- extract Proba of both options
        probaJ1( when('I_juiceL',1) & restr ) = TrialType(take('I_probaL'),when('I_juiceL',1) & restr);
        probaJ1( when('I_juiceR',1) & restr ) = TrialType(take('I_probaR'),when('I_juiceR',1) & restr);
        probaJ1 = probaJ1(restr);
        probaJ2( when('I_juiceL',2) & restr ) = TrialType(take('I_probaL'),when('I_juiceL',2) & restr);
        probaJ2( when('I_juiceR',2) & restr ) = TrialType(take('I_probaR'),when('I_juiceR',2) & restr);
        probaJ2 = probaJ2(restr);

        pbJ1 = norm(probaJ1/100); %- normalize proba
        pbJ2 = norm(probaJ2/100);

        %- extract Reaction time
        rt = t_respfix-t_resp(1,:);
        rt = norm_RT(log(rt(restr))); %- log then min-max normalization (-1 to 1)
        rt(~isfinite(rt)) = NaN;

        %- extract Initiation time
        ft = t_FPfix_first-t_FPon(1,:);
        ft = norm_FT(log(ft(restr))); %- log then min-max normalization (-1 to 1)
        ft(~isfinite(ft)) = NaN;

        %- extract prevJuice (0 if chose J1 before / 1 if chose J2 before) - this is not side related
        %- prevJuice == which juice he 'could have get' on the last completed trial (not necessarily
        %- rewarded).
        %- WARNING!!! FOR NOW, INCLUDE Dynamic task trials -- If trial before is Dynamic task, which juice he chose there (Morbier) or was offered (Mimic)

        currjuice = NaN(1,length(TrialType));
        currjuice(when('task',1)) = TrialType(take('P_juice'),when('task',1))-1;
        currjuice(when('task',2)) = TrialType(take('I_chosenjuice'),when('task',2))-1;
        currjuice(when('task',5)) = TrialType(take('D_chosenjuice'),when('task',5))-1;
        currjuice(when('task',3)) = TrialType(take('D_chosenjuice'),when('task',3))-1;

        %- put the last trials' juice when brks at current (or multiple breaks in a row, take the last one)
        noJuice = find(isnan(currjuice));
        for tr = 1 : length(noJuice)
            if noJuice(tr)>1
                currjuice(noJuice(tr)) = currjuice(noJuice(tr)-1);
            end
        end

        prevJuice = [NaN currjuice(1:end-1)];
        prevJuice = prevJuice(restr);
        prevJuice(prevJuice==0)=-1; %- a way to standardize it (-1 vs 1)

        %- extract the last received juice
        received = when('rew',1);
        currjuice(~received)=NaN; %- remove all the juice info from unrewarded trials

        %- put the last trials' RECEIVED juice when no rew or brks (or multiple breaks in a row, take the last one)
        %- no information on how far was this last reward
        noRwd = find(isnan(currjuice));
        for tr = 1 : length(noRwd)
            if noRwd(tr)>1
                currjuice(noRwd(tr)) = currjuice(noRwd(tr)-1);
            end
        end

        prevJuiceReceived = [NaN currjuice(1:end-1)];
        prevJuiceReceived = prevJuiceReceived(restr);
        prevJuiceReceived(prevJuiceReceived==0)=-1; %- a way to standardize it (-1 vs 1)

        %- define model
        % modelspec = 'choice ~ 1 + probJ1*probJ2 + ft + prevJuice' ;
        modelspec = 'choice ~ 1 + probJ1*probJ2 + prevJuice' ;
        % modelspec = 'choice ~ 1 + probJ1*probJ2 + prevJuiceReceived' ;
        % modelspec = 'choice ~ 1 + probJ1*probJ2 + rt' ;

        %- create Table and set categorical variables
        T = table(choice',pbJ1',pbJ2',ft',rt',prevJuice',prevJuiceReceived','VariableNames',{'choice','probJ1','probJ2','ft','rt','prevJuice','prevJuiceReceived'});
        T.choice=categorical(T.choice);
        T.prevJuice=categorical(T.prevJuice);
        T.prevJuiceReceived=categorical(T.prevJuiceReceived);

        lastwarn('', '');
        mdl = fitglm(T,modelspec,'Distribution','binomial','Link','logit')

        %- check if glm converged
        [~, warnId] = lastwarn();
        if strcmp(warnId,'stats:glmfit:IterationLimit')
            converge=0;
        else
            converge=1;
        end

        %- can't really test for normality given the number of samples.. any test like KS will be
        %- significant with large sample size.. also see : web.pdx.edu/~newsomj/cdaclass/ho_diagnostics.pdf
        % figure;qqplot(mdl.Residuals.Raw)
        % figure;normplot(mdl.Residuals.Raw)
        % [h,p] = kstest(mdl.Residuals.Raw)

        %- extract residuals
        residuals = mdl.Residuals.Raw;
        residuals_NaNfree = residuals;
        residuals_NaNfree(isnan(residuals_NaNfree))=[];

        %- test for trend in residuals // monotonic trend, but not necessarily linear
        [Z_trend,p_trend]=Mann_Kendall(residuals_NaNfree,0.01)

        %- predict choices (for that, prevJuice is set as 0 (J1) and RT is median(rt)
        clear newf newc
        probas = norm([10:1:90]);
        thr = NaN(length(probas),1);

        for i = 1 : length(probas)

            tab = table((probas(i)*ones(length(probas),1)),probas',nanmedian(ft)*ones(length(probas),1),nanmedian(rt)*ones(length(probas),1),-ones(length(probas),1),-ones(length(probas),1),'VariableNames',{'probJ1','probJ2','ft','rt','prevJuice','prevJuiceReceived'});
            tab.prevJuice = categorical(tab.prevJuice);
            tab.prevJuiceReceived = categorical(tab.prevJuiceReceived);

            [newf(i,:) , newc] = predict(mdl,tab); %- last column is the sessiooon

            if ~isempty(find(newf(i,:)>=0.5,1,'first'))
                thr(i,1) = find(newf(i,:)>=0.5,1,'first');
            end
        end

        %- extract the preference for J1 vs J2
        %- sensitivity doesn't take into account the interaction factor... so can use the nb of
        %- estimated probaChoiceJ1>0.5.. if pref>0.5, prefer the J1
        sensitivity = mdl.Coefficients{'probJ1','Estimate'}/mdl.Coefficients{'probJ2','Estimate'};
        pref = sum(sum(newf<0.5))/(length(newf(:))-sum(sum(newf==0.5)));

        %- confidence (can use that as regressor of neuronal activity
        conf=-newf.*log(newf)-(1-newf).*log(1-newf); %% define by amemori & graybiel 2012 last supp figure


        clear sensitivity_bins pref_bins radj_bins ft_bins rt_bins mdl_bins resid
        %     bins = floor(height(T)/5);
        %
        %         bins_start = 1:floor(height(T)-bins/12):height(T)-bins;
        % add some -1?!?!? doesn't work now...
        bins = floor(height(T)/5);
        bins_start = round(1:(height(T)-bins)/50:height(T)-bins);
        parfor b = 1 : length(bins_start)
            warning off
            disp(b)
            mdl_bins = fitglm(T(bins_start(b):bins_start(b)+bins,:),modelspec,'Distribution','binomial','Link','logit');
            ft_bins(b) = median(table2array(T(bins_start(b):bins_start(b)+bins,'ft')));
            rt_bins(b) = median(table2array(T(bins_start(b):bins_start(b)+bins,'rt')));
            resid(1,b) = mean(mdl.Residuals.Raw(bins_start(b):bins_start(b)+bins));

            for i = 1 : length(probas)
                tab = table(probas(i)*ones(length(probas),1),probas',nanmedian(ft)*ones(length(probas),1),nanmedian(rt)*ones(length(probas),1),-ones(length(probas),1),-ones(length(probas),1),'VariableNames',{'probJ1','probJ2','ft','rt','prevJuice','prevJuiceReceived'});
                tab.prevJuice = categorical(tab.prevJuice);
                tab.prevJuiceReceived = categorical(tab.prevJuiceReceived);
                if i == 1
                    tab_all = tab;
                else
                    tab_all = vertcat(tab_all,tab);
                end
            end
            newf_bins = [];
            try
                [newf_bins , ~] = predict(mdl_bins,tab_all);
                newf_bins = reshape(newf_bins,length(probas),length(probas));
            end

            sensitivity_bins(b) = mdl_bins.Coefficients{'probJ1','Estimate'}/mdl_bins.Coefficients{'probJ2','Estimate'};
            if ~isempty(newf_bins)
                pref_bins(b) = sum(sum(newf_bins<0.5))/(length(newf_bins(:))-sum(sum(newf_bins==0.5)));
            else
                pref_bins(b) = NaN;
            end
            radj_bins(b) = mdl_bins.Rsquared.Adjusted;
        end
        %- test for trend in pref bins
        [Z_trend_bins,p_trend_bins]=Mann_Kendall(pref_bins,0.01);
        [r_corr,p_corr] = corrcoef(ft_bins,pref_bins);

        %% Now check behavior on same juice trials

        %- check the choice consistency in same juice trials (in individually for each Juices)
        %- choice consistency is how often he selected the higher proba!
        same_juice = (when('I_juiceL',1) & when('I_juiceR',1)) | (when('I_juiceL',2) & when('I_juiceR',2));
        restr_same = when('task',2) & when('brks',0) & same_juice; %- take only INS trials (task == 2), completed (brks == 0) and same juice

        chosenPb_J1 = TrialType(take('I_chosenproba'),restr_same & when('I_juiceL',1));
        unchosenPb_J1 = TrialType(take('I_unchosenproba'),restr_same & when('I_juiceL',1));
        chosenPb_J2 = TrialType(take('I_chosenproba'),restr_same & when('I_juiceL',2));
        unchosenPb_J2 = TrialType(take('I_unchosenproba'),restr_same & when('I_juiceL',2));
        chosenPb = TrialType(take('I_chosenproba'),restr_same);
        unchosenPb = TrialType(take('I_unchosenproba'),restr_same);

        choice_consistency_J1 = sum(chosenPb_J1>unchosenPb_J1)/length(chosenPb_J1);
        choice_consistency_J2 = sum(chosenPb_J2>unchosenPb_J2)/length(chosenPb_J2);

        %- check RTs in same juice trials
        rt_same = t_respfix-t_resp(1,:);
        rt_same = norm(log(rt_same(restr_same))); %- min-max normalization + log
        rt_same(~isfinite(rt_same)) = NaN;

        %- check RTs in same juice trials
        ft_same = t_FPfix_first-t_FPon(1,:);
        ft_same = norm(log(ft_same(restr_same))); %- min-max normalization + log
        ft_same(~isfinite(ft_same)) = NaN;

        %- compute choice consistency on the same bins from the model before
        nTr = 1:length(TrialType);
        diffIns_tr = nTr(restr);
        sameIns_tr = nTr(restr_same);
        clear choice_consistency_bins ft_same_bins rt_same_bins
        for b = 1 : length(bins_start)
            boundaries = diffIns_tr([bins_start(b) bins_start(b)+bins]);
            tr_bins = sameIns_tr>=boundaries(1) & sameIns_tr<=boundaries(2);
            choice_consistency_bins(b) = sum(chosenPb(tr_bins)>unchosenPb(tr_bins))/length(chosenPb(tr_bins));
            ft_same_bins(b) = median(ft_same(tr_bins));
            rt_same_bins(b) = median(rt_same(tr_bins));
        end

        %%
        mat_name = {'Residuals' 'Juice pref' 'IT' 'RT' 'IT same' 'RT same' 'Consistency same'};
        mat = [resid ; pref_bins ; ft_bins ; rt_bins ; ft_same_bins ; rt_same_bins ; choice_consistency_bins];
        for p1 = 1 : length(mat(:,1))
            for p2 = 1 : length(mat(:,1))
                if p1>p2
                    [r,p] = corrcoef(mat(p1,:),mat(p2,:));
                    r_all(p1,p2) = r(2,1);
                    p_all(p1,p2) = p(2,1);
                else r_all(p1,p2) = NaN;
                    p_all(p1,p2) = NaN;
                end
            end
        end

        %% Extract param of interest for each session!
        ALL(dd,1).name = filename;
        ALL(dd,1).mdl = mdl;
        ALL(dd,1).converge = converge;
        ALL(dd,1).resid = resid;
        ALL(dd,1).pref = pref;
        ALL(dd,1).Z_trend = Z_trend;
        ALL(dd,1).p_trend = p_trend;
        ALL(dd,1).p_trend_bins = p_trend_bins;
        ALL(dd,1).Z_trend_bins = Z_trend_bins;
        ALL(dd,1).pref_bins = single(pref_bins);
        ALL(dd,1).radj_bins = single(radj_bins);
        ALL(dd,1).ft_bins = single(ft_bins);
        ALL(dd,1).rt_bins = single(rt_bins);
        ALL(dd,1).choice_consistency_bins = single(choice_consistency_bins);
        ALL(dd,1).ft_same_bins = single(ft_same_bins);
        ALL(dd,1).rt_same_bins = single(rt_same_bins);
        ALL(dd,1).corrmat = r_all;
        ALL(dd,1).corrmat_p = p_all;

        if ~perm %- keep more param when not perm
            ALL(dd,1).TrialType = single(TrialType);
            ALL(dd,1).trial2take = restr;
            ALL(dd,1).T = T;
            ALL(dd,1).sensitivity = sensitivity;
            ALL(dd,1).conf = single(conf);
            ALL(dd,1).rt = single(rt);
            ALL(dd,1).ft = single(ft);
            ALL(dd,1).r_corr = r_corr;
            ALL(dd,1).p_corr = p_corr;
            ALL(dd,1).sensitivity_bins = single(sensitivity_bins);
            ALL(dd,1).choice_consistency_J1 = choice_consistency_J1;
            ALL(dd,1).choice_consistency_J2 = choice_consistency_J2;
            ALL(dd,1).rt_same = single(rt_same);
            ALL(dd,1).ft_same = single(ft_same);
        end
        if perm
            ALL(dd,1).tr_order = tr_order; %- perm the others (ins + dyn)
        end

        param.modelspec = modelspec;
        param.binSize = bins;
        param.mat_name = mat_name;

        %- plot a quick recap figure, if needed

        if showfig

            figure;
            colors = cbrewer('div', 'PiYG', 64);
            colors = flipud(colors); % puts pink on top, green at the bottom

            subplot(3,4,1);
            imagesc(newf,[0 1]);axis xy
            colormap(colors);
            hold on
            plot(thr(:,1),1:length(probas),'-k','LineWidth',2);ylabel('Proba J1');xlabel('Proba J2');

            subplot(3,4,2)
            imagesc(conf);axis xy
            colormap(colors);
            title([mdl.Formula.LinearPredictor ' / BIC: ' num2str(mdl.ModelCriterion.BIC) ' / Disp= ' num2str(mdl.Dispersion) ' / Pref= ' num2str(pref) ' / R2= ' num2str(mdl.Rsquared.Adjusted)  ])

            subplot(3,4,[3 4]);
            bar(mdl.Coefficients.Estimate,'FaceColor',[.6 .6 .6])
            set(gca,'XTick',1:length(mdl.Coefficients.Estimate),'XTickLabels',mdl.CoefficientNames,'XTickLabelRotation',20	)
            pval =  round(mdl.Coefficients.pValue*10000)/10000;
            for i = 1 : length(mdl.Coefficients.tStat)
                if pval(i)<0.01
                    text(i,mdl.Coefficients.Estimate(i),num2str(pval(i)),'Color','r','HorizontalAlignment','center')
                else
                    text(i,mdl.Coefficients.Estimate(i),num2str(pval(i)),'HorizontalAlignment','center')
                end
            end
            xlim([0.5 length(pval)+.5])

            subplot(3,4,[5 6]);
            if p_trend<0.01 ; plot(mdl.Residuals.Raw,'r');
            else plot(mdl.Residuals.Raw,'k');
            end
            text(25,-0.75,['p-trend=' num2str(p_trend)])
            ylabel('Residuals');xlabel('Trials');

            subplot(3,4,[9 10]);
            line([0 length(pref_bins)],[0.5 0.5],'Color','k');box on;hold on
            if p_trend_bins<0.01 ; plot(pref_bins,'r')
            else plot(pref_bins,'k')
            end
            text(25,0.15,['p-trend=' num2str(p_trend_bins)])
            ylabel('Preference (J1 when > 0.5)');
            hold off
            ylim([0 1])

            subplot(3,4,[11]);
            if p_corr(1,2)<0.01
                plot(ft_bins,pref_bins,'.r')
            else
                plot(ft_bins,pref_bins,'.k')
            end
            xlabel('FT')
            title(['R=' num2str(r_corr(1,2)) ' / p=' num2str(p_corr(1,2))])
            ylim([0 1])

        end

    else
        ALL(dd,1) = [];
        param(dd,1) = [];
    end

end