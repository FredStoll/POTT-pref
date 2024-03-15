%% POTT Pref - Decoding (LDA) on pseudo-population of neurons (Figure S1/S2)
%-
%- Works for running individual monkeys and combined, for both full set and
%- subset of neurons.
%-
%- if you want to run additional predictors, it will update the previously
%- generated matrix.
%-
%- Author: Fred M. Stoll, Icahn School of Medicine at Mount Sinai, NY
%- Date: 2023.05
%- Related to: Stoll & Rudebeck, Neuron, 2024

clear
monkey2run = {'M' 'X' 'B'}; % 'M', 'X', 'Both' ('Both' = combined)
subsetNeurons = true;
% path2go = '/home/fred/POTT/data/subset-final/'; %- path where SPKpool files are!
path2go = '/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-BehavPrefChange/data/neurons/subset-final/'; %- path where SPKpool files are!

for mmm = 1 : length(monkey2run)

    clearvars -except monkey2run mmm path2go subsetNeurons

    %- define which file to use (subset of neurons or full, and # of neurons to run)
    if subsetNeurons
        startfilename = 'POOLsub';
        savename = 'res_LDApopsub_kfold_final';
        param.minNeurons = [50 100 150 200 250 300 400 500]; %- min number of neuron to run
    else
        startfilename = 'POOL';
        savename = 'res_LDApop_kfold_final';
        param.minNeurons = [50 100 150 200 250 300 400 500 750 1000 1250 1500 2000 2500]; %- min number of neuron to run
    end

    %- locate files
    if strcmp(monkey2run{mmm},'M')
        list = dir([path2go startfilename '_97*.mat']);
    elseif strcmp(monkey2run{mmm},'X')
        list = dir([path2go startfilename '_186*.mat']);
    elseif strcmp(monkey2run{mmm}(1),'B')
        list = dir([path2go startfilename '_283*.mat']);
    end
    params = load([path2go 'M021519a_SPKpool.mat']); %- to get the param from a random recording session

    % all the param required
    param.predic = {'I_chosenjuice' 'I_chosenproba' 'I_chosenside' 'I_unchosenproba'}; %- select param you want to test..
    param.nComp = 20; %- nb of component to keep for PCA
    param.Repetition = 200; %- number of times to run the decoder (100 or 200 ideally)
    param.minTrCond = 10; %- minimum number of trials per condition to run
    param.overwrite = false;
    param.thresh = false;
    param.thr = 1; % in Hz, used if param.thresh is 'y'
    param.rmv = 100; % remove xxx ms on each side on every events (avoid overlaps between bins and smoothing problems)
    param.bins4decoding = 2 ; %- 6 for rew %- perform decoding on subset on bins (stim period here)
    param.timebin = [200 800] ;% 0 600 rew %- perform decoding on subset on bins (stim period here)
    param.normalize_fr = false;
    param.tr_perm = true;

    disp(['Computing LDA on ' num2str(length(list)) ' sessions'])
    areas = utils_POTT_areas;

    area2test = {'vlPFC' 'OFC' 'IFG' 'LAI' 'AMG' };

    util.minmaxnorm = @(data) (data-min(data))/(max(data)-min(data));
    util.meannorm  = @(data) (data-mean(data))/(max(data)-min(data));
    util.minnorm  = @(data) (data-min(min(data)))/(max(max(data))-min(min(data)));

    un = find(list(1).name=='_');

    %- check if matrix already exists for each predictor
    if exist([path2go savename '_' list(1).name(un(1)+1:un(2)-1) '.mat'])==2
        prev = load([path2go savename '_' list(1).name(un(1)+1:un(2)-1) '.mat']);
        done = [];
        for pr = 1 : length(param.predic)
            if isfield(prev.res_LDApop,param.predic{pr})
                done(pr) = true;
            end
        end
        if ~param.overwrite
            param.predic(find(done==true))=[];
        end
        res_LDApop = prev.res_LDApop;
    end
    if param.overwrite
        clear res_LDApop
    end

    if ~isempty(param.predic) %- if there is some predictors not done yet

        [~,name] = system('hostname');
        if isempty(gcp('nocreate')) && strcmp(name(1:5),'DESKT')
            parpool(8);
        elseif isempty(gcp('nocreate')) && ~strcmp(name(1:5),'DESKT')
            parpool(36);
        end

        for pr = 1 : length(param.predic)
            %nSess = zeros(length(area2test),1);
            for s = 1 : length(list)
                clearvars -except path2go list param params areas area2test util s u pr res_LDApop nSess un
                name = list(s).name;

                %- load spike data and histology info
                load([path2go name]); disp(s);


                %- create time vector and remove overlapping windows/border events
                bins2remove = (param.rmv/params.subsp)-1;

                n_evt = length(params.times_evts);
                lin=size(neurons_info,1)*n_evt;
                nTimes = [0 ; (sum(abs(params.times(:,:)),2)/params.subsp)];
                col = sum(nTimes);
                bins = NaN(1,col);
                time = NaN(1,col);

                for b = 1 : length(params.times_evts)
                    time(sum(nTimes(1:b))+1:sum(nTimes(1:b+1))) =  (params.times(b,1):params.subsp:params.times(b,2)-(params.subsp/2));
                    bins(sum(nTimes(1:b))+1:sum(nTimes(1:b+1))) = b;
                    bins([sum(nTimes(1:b))+1:sum(nTimes(1:b))+1+bins2remove   sum(nTimes(1:b+1))-bins2remove:sum(nTimes(1:b+1))    ]) = NaN; %- remove 100 ms each side (avoid overlaps between bins and smoothing problems)
                end

                %- extract the data for the task 2 consider
                if strcmp(param.predic{pr}(1),'I') %- instrumental
                    SPK_data = SPK_INS;
                    BEHAV = BEHAV_INS;
                end

                %- normalize FR + cut to only bin 4 decoding
                if param.normalize_fr %- normalize FR using min-max in considered bins only (all trials except border of events)
                    SPK_data_norm = SPK_data;
                    units = unique(Tr_Clust_data(:,2));
                    for u = 1 : length(units)
                        temp = SPK_data(Tr_Clust_data(:,2)==units(u),~isnan(bins)); %- could normalize only on the bin2use, but more likely to get some 0 FR?? so NaNs!
                        SPK_data_norm(Tr_Clust_data(:,2)==units(u),~isnan(bins)) = reshape(util.minnorm(temp(:)),size(temp,1),size(temp,2));
                    end

                    SPK_data_cut = SPK_data_norm;
                    SPK_data_cut(:,~ismember(bins,param.bins4decoding))=[];
                else
                    SPK_data_cut = SPK_data;
                    SPK_data_cut(:,~ismember(bins,param.bins4decoding))=[];
                end

                time(~ismember(bins,param.bins4decoding)) = [];

                %- reject low FR neurons, if needed
                if param.thresh  %- reject neurons with FR too low
                    reject = [];
                    for n = 1 : length(neurons_rank)
                        if mean(mean(SPK_data_cut(Tr_Clust_data(:,2)==n,:)))<param.thr  %- put NaNs when FR is below thr
                            reject = [reject , n];
                        end
                    end

                    %- remove them from main matrices
                    SPK_data_cut( ismember(Tr_Clust_data(:,2),reject) ,:) = [];
                    Tr_Clust_data( ismember(Tr_Clust_data(:,2),reject) ,:) = [];
                    area_histology(reject)=[];
                end

                %- create data and factor matrix for the decoding
                data = SPK_data_cut;
                unit_nb = unique(BEHAV(1,:));

                TrialType_header = ['Unit_Nb' 'Session_Nb' params.TrialType_header];

                %- trials to consider for instrum
                if strcmp(param.predic{pr}(1),'I')
                    when = @(arg,cd) BEHAV(strcmp(TrialType_header,arg),:)==cd ;
                    diff_juice = (when('I_juiceL',1) & when('I_juiceR',2)) | (when('I_juiceL',2) & when('I_juiceR',1)); %- take only the different juice trials
                end

                if strcmp(param.predic{pr}(1:end-1),'I_offerproba_j')
                    offerJ1 = NaN(size(diff_juice));
                    offerJ2 = NaN(size(diff_juice));

                    offerJ1(when('I_juiceL',1) & diff_juice) = BEHAV(ismember(TrialType_header,'I_probaL'),when('I_juiceL',1) & diff_juice) ;
                    offerJ1(when('I_juiceR',1) & diff_juice) = BEHAV(ismember(TrialType_header,'I_probaR'),when('I_juiceR',1) & diff_juice );
                    offerJ2(when('I_juiceL',2) & diff_juice) = BEHAV(ismember(TrialType_header,'I_probaL'),when('I_juiceL',2) & diff_juice) ;
                    offerJ2(when('I_juiceR',2) & diff_juice) = BEHAV(ismember(TrialType_header,'I_probaR'),when('I_juiceR',2) & diff_juice );

                    dumm = BEHAV(ismember(TrialType_header,{'Unit_Nb' }),:)';
                    eval(['factor = [transpose(offerJ' param.predic{pr}(end) ') dumm];'])
                    dumm = factor(:,[2 1]);
                else
                    dumm = BEHAV(ismember(TrialType_header,{'Unit_Nb' param.predic{pr} }),:)';
                    factor = dumm(:,[2 1]);
                end

                %- only keep a subset of proba
                remove = [];
                if strcmp(param.predic{pr},'P_proba')
                    keep = [10 30 50 70 90 ];
                    remove = ~ismember(dumm(:,2),keep);
                elseif strcmp(param.predic{pr},'I_chosenproba')
                    keep = [30 50 70 90 ];
                    remove = ~ismember(dumm(:,2),keep) | ~diff_juice';  % was not removing same juice trials before
                elseif strcmp(param.predic{pr},'I_unchosenproba')
                    keep = [10 30 50 70];
                    remove = ~ismember(dumm(:,2),keep) | ~diff_juice';  % was not removing same juice trials before
                elseif strcmp(param.predic{pr},'I_chosenjuice')
                    keep = [1 2];
                    remove = ~diff_juice';  % was not removing same juice trials before
                elseif strcmp(param.predic{pr}(1:end-1),'I_offerproba_j')
                    keep = [10 30 50 70 90 ];
                    remove = ~ismember(dumm(:,2),keep) | ~diff_juice';  % was not removing same juice trials before
                else
                    keep = unique(dumm(:,2));
                end
                if ~isempty(remove)
                    data(remove,:) = [];
                    dumm(remove,:)=[];
                    factor(remove,:)=[];
                end

                %- check if enough trial and remove neurons without enough
                for p = 1 : length(keep)
                    for n = 1 : length(unit_nb)
                        nTr(n,p) = sum(dumm(:,2)==keep(p) & dumm(:,1)==unit_nb(n));
                    end
                end

                remove = sum(nTr>param.minTrCond,2)~= length(keep);

                unit_nb(remove)=[];

                %- extract only 1 bin during stim onset
                data = mean(data(:,time>=param.timebin(1) & time<=param.timebin(2)),2);

                for u = 1 : length(param.minNeurons)

                    if length(unit_nb)>=param.minNeurons(u)

                        parfor p = 1 : param.Repetition
                            disp(['Area ' num2str(s) '/' num2str(length(list)) ' - ' num2str(param.minNeurons(u)) ' units - ' param.predic{pr} ' - perm ' num2str(p) '/' num2str(param.Repetition)])

                            if ~isempty(param.minNeurons(u))
                                if param.minNeurons(u)<=length(unit_nb)
                                    rd = randperm(length(unit_nb));
                                    unit_nb_sub = unit_nb(rd);
                                    unit_nb_sub = unit_nb_sub(1:param.minNeurons(u));
                                end
                            end

                            data_sub =  data( ismember(factor(:,2),unit_nb_sub),:);
                            factor_sub =  factor( ismember(factor(:,2),unit_nb_sub),:);

                            [XX,YY,param_LDA] = pop_init_noRdmTrials(data_sub,factor_sub,1,'perm',1,'minTr',[],'pop',false); % ,'window',[-.5 1.5]
                            [perf(:,p),~,out_dumm] = pca_lda_kfold(XX,YY,param);
                            out(p).recall = out_dumm.recall;
                            out(p).precision=out_dumm.precision;
                            out(p).specificity=out_dumm.specificity;
                            out(p).f_meas=out_dumm.f_meas  ;

                            %- permute trial factors
                            if param.tr_perm
                                loopi = unique(factor_sub(:,2));
                                for up = 1 : length(loopi)
                                    dumm_condi = factor_sub(factor_sub(:,2)==loopi(up),1);

                                    factor_sub(factor_sub(:,2)==loopi(up),1) = dumm_condi(randperm(length(dumm_condi)));
                                end

                                [XXperm,YYperm,param_LDA] = pop_init_noRdmTrials(data_sub,factor_sub,1,'perm',1,'minTr',[],'pop',false); % ,'window',[-.5 1.5]
                                [perf_perm(:,p),~,out_perm_dumm] = pca_lda_kfold(XXperm,YYperm,param);
                                out_perm(p).recall = out_perm_dumm.recall;
                                out_perm(p).precision=out_perm_dumm.precision;
                                out_perm(p).specificity=out_perm_dumm.specificity;
                                out_perm(p).f_meas=out_perm_dumm.f_meas  ;
                            end


                        end

                        %- extract param of interest to save
                        %  nSess(ar,1) = nSess(ar,1) + 1 ;
                        temp.lda_sess = name;
                        temp.perf = perf;
                        temp.nbUnit = param.minNeurons(u);
                        temp.output = out;
                        if param.tr_perm
                            temp.perf_perm = perf_perm;
                            temp.output_perm = out_perm;
                        end
                        unde = find(name=='_',1,'last');
                        area_name = ['area_' name(unde+1:end-4)];
                        eval(['res_LDApop.' param.predic{pr} '.' area_name '(u) = temp;'])
                    end
                end
            end
            eval(['res_LDApop.' param.predic{pr} '.param = param;'])
        end

        save([path2go savename '_' list(1).name(un(1)+1:un(2)-1) '.mat'],'res_LDApop','list','area2test','path2go')
    end

end

%% POST HOC - COMBINED MONKEY - FULL DATASET (Fig S1)

clear
path2go = 'S:\POTT\data\subset-final\'; %- path where SPKpool files are!
load([path2go 'res_LDApop_kfold_final_283sessions.mat'])
measures = {'I_chosenjuice' 'I_chosenproba' 'I_chosenside'};
nbNeur = 500;

[colorsArea,colorsArea_sub] = colorMeUp('POTT'); %- COLOR ASSIGNMENT

clear recall precision specificity f_meas
clear recall_perm precision_perm specificity_perm f_meas_perm
meanPerf={};
meanPerf_perm={};
nbNeur_sig =[];
for m = 1 : length(measures)
    x=0;allPerf{m}=[];
    for ar = 1: length(area2test)
        eval(['temp = res_LDApop.' measures{m} '.area_' area2test{ar} ';'])
        for n = 1 : length(temp)
            %      if temp(n).nbUnit<=500
            x = x+1;
            pval = sum(temp(n).perf_perm > nanmean(temp(n).perf)) / length(temp(n).perf_perm);
            sem = (nanstd(temp(n).perf)/sqrt(length(temp(n).perf)));

            recall{m}(x,:) = nanmean(cat(1,temp(n).output(:).recall)   );
            precision{m}(x,:) = nanmean(cat(1,temp(n).output(:).precision)   );
            specificity{m}(x,:) = nanmean(cat(1,temp(n).output(:).specificity)   );
            f_meas{m}(x,:) = nanmean(cat(1,temp(n).output(:).f_meas)   );

            recall_perm{m}(x,:) = nanmean(cat(1,temp(n).output_perm(:).recall)   );
            precision_perm{m}(x,:) = nanmean(cat(1,temp(n).output_perm(:).precision)   );
            specificity_perm{m}(x,:) = nanmean(cat(1,temp(n).output_perm(:).specificity)   );
            f_meas_perm{m}(x,:) = nanmean(cat(1,temp(n).output_perm(:).f_meas)   );

            meanPerf{m}(x,:) = [nanmean(temp(n).perf) temp(n).nbUnit ar nanmean(temp(n).perf)-sem nanmean(temp(n).perf)+sem pval];
            sem_perm = (nanstd(temp(n).perf_perm)/sqrt(length(temp(n).perf_perm)));
            meanPerf_perm{m}(x,:) = [nanmean(temp(n).perf_perm) temp(n).nbUnit ar nanmean(temp(n).perf_perm)-sem_perm nanmean(temp(n).perf_perm)+sem_perm];
            nP = length(temp(n).perf);
            allPerf{m} = [allPerf{m} ; temp(n).perf' repmat(temp(n).nbUnit(1),nP,1) repmat(ar,nP,1) (1:nP)'];
        end
        % end
    end

    for ar = 1: length(area2test)
        idx = find(meanPerf{m}(:,3)==ar & meanPerf{m}(:,end)<0.05,1,'first');
        if ~isempty(idx)
            nbNeur_sig(m,ar) = meanPerf{m}(idx,2);
        else
            nbNeur_sig(m,ar) = NaN;
        end
        idx2 = find(meanPerf{m}(:,3)==ar ,1,'last');
        nbNeur_test(m,ar) = meanPerf{m}(idx2,2);
    end
end

probas = {'30%' '50%' '70%' '90%'};
figure;m = 2;
subplot(2,2,1);
for ar = 1: length(area2test)
    dumm = [recall{m}(meanPerf{m}(:,3)==ar & meanPerf{m}(:,2)==nbNeur,:)];
    plot(dumm,'.-','Color',colorsArea(ar,:),'LineWidth',2,'MarkerSize',15);hold on
    dumm = [recall_perm{m}(meanPerf{m}(:,3)==ar & meanPerf{m}(:,2)==nbNeur,:)];
    plot(dumm,'--','Color',colorsArea(ar,:),'LineWidth',2,'MarkerSize',15);hold on
end
set(gca,"XTick",1:4,'XTickLabel',probas);xlim([0 5])
ylim([0 1]);title('Recall')
subplot(2,2,2);
for ar = 1: length(area2test)
    dumm = [precision{m}(meanPerf{m}(:,3)==ar & meanPerf{m}(:,2)==nbNeur,:)];
    plot(dumm,'.-','Color',colorsArea(ar,:),'LineWidth',2,'MarkerSize',15);hold on
    dumm = [precision_perm{m}(meanPerf{m}(:,3)==ar & meanPerf{m}(:,2)==nbNeur,:)];
    plot(dumm,'--','Color',colorsArea(ar,:),'LineWidth',2,'MarkerSize',15);hold on
end
set(gca,"XTick",1:4,'XTickLabel',probas);xlim([0 5])
ylim([0 1]);title('Precision')
subplot(2,2,3);
for ar = 1: length(area2test)
    dumm = [specificity{m}(meanPerf{m}(:,3)==ar & meanPerf{m}(:,2)==nbNeur,:) ];
    plot(dumm,'.-','Color',colorsArea(ar,:),'LineWidth',2,'MarkerSize',15);hold on
    dumm = [specificity_perm{m}(meanPerf{m}(:,3)==ar & meanPerf{m}(:,2)==nbNeur,:) ];
    plot(dumm,'--','Color',colorsArea(ar,:),'LineWidth',2,'MarkerSize',15);hold on
end
set(gca,"XTick",1:4,'XTickLabel',probas);xlim([0 5])
ylim([0.5 1]);title('Specificity')
subplot(2,2,4);
for ar = 1: length(area2test)
    dumm = [f_meas{m}(meanPerf{m}(:,3)==ar & meanPerf{m}(:,2)==nbNeur,:) ];
    plot(dumm,'.-','Color',colorsArea(ar,:),'LineWidth',2,'MarkerSize',15);hold on
    dumm = [f_meas_perm{m}(meanPerf{m}(:,3)==ar & meanPerf{m}(:,2)==nbNeur,:) ];
    plot(dumm,'--','Color',colorsArea(ar,:),'LineWidth',2,'MarkerSize',15);hold on
end
set(gca,"XTick",1:4,'XTickLabel',probas);xlim([0 5])
ylim([0 1]);title('F-measure (~ Recall/Precision ratio)')

%- stats of info and total infor
chance_l = [0.5 0.25 0.5];

for m = 1 : length(measures)
    for ar = 1: length(area2test)
        nPerm = length(unique(allPerf{m}(:,4)));
        for p = 1 : nPerm
            dumm = allPerf{m}(allPerf{m}(:,3)==ar & allPerf{m}(:,4)==p,[1 2]);
            [param_a{m}(ar,p),param_b{m}(ar,p),~,exitflag{m}(ar,p)] = fitDecodPerf(dumm,chance_l(m));
        end
    end
    param_b{m} = param_b{m}+chance_l(m);
    param_b{m}(exitflag{m}==0)=NaN;
    param_a{m}(exitflag{m}==0)=NaN;
end

for m = 1 : length(measures)
    dumm_areas=repmat((1:length(area2test))',1,length(param_b{m}));
    modeldata_all = table(param_a{m}(:),param_b{m}(:),area2test(dumm_areas(:))','VariableNames',{'a','b','area'});

    models_form = {'a ~ 1 + area'};
    lme_a = fitglme(modeldata_all(~isnan(modeldata_all.a),:),models_form{1});
    [pval_a,wald_a,thr_corr,pval_a_adj] = area_posthoc(lme_a,area2test,'y');
    disp(anova(lme_a))
    models_form = {'b ~ 1 + area'};
    lme_b = fitglme(modeldata_all(~isnan(modeldata_all.b),:),models_form{1});
    [pval_a,wald_a,thr_corr,pval_a_adj] = area_posthoc(lme_b,area2test,'y');
    disp(anova(lme_b))
end

%- MAIN FIGURE
subplot_idx = {[ 1 2  7  8] 13 14 ;
    [ 3 4  9 10] 15 16 ;
    [ 5 6 11 12] 17 18 };
chance_l = [0.5 0.25 0.5];

figure;
for m = 1 : length(measures)
    subplot(3,6,subplot_idx{m,1});
    for ar = 1: length(area2test)
        dumm = meanPerf{m}(meanPerf{m}(:,3)==ar,[1 2 4 5]);
        plot(dumm(:,2),dumm(:,1),'.-','Color',colorsArea(ar,:),'LineWidth',2,'MarkerSize',15);hold on

        [~,~,Yh] = fitDecodPerf(dumm,chance_l(m));
        plot([0 ; dumm(:,2)],[0 ; Yh]+chance_l(m),'--','Color',colorsArea(ar,:),'LineWidth',1,'MarkerSize',15);hold on

        for p = 1:size(dumm,1)
            line([dumm(p,2) dumm(p,2)],[dumm(p,3) dumm(p,4)],'Color',colorsArea(ar,:))
        end
        keep4legend(ar,1) = dumm(1,1);

        dumm = meanPerf_perm{m}(meanPerf_perm{m}(:,3)==ar,[1 2 4 5]);
        plot(dumm(:,2),dumm(:,1),'.-','Color',colorsArea(ar,:),'LineWidth',1,'MarkerSize',10);hold on
        for p = 1:size(dumm,1)
            line([dumm(p,2) dumm(p,2)],[dumm(p,3) dumm(p,4)],'Color',colorsArea(ar,:))
        end

        % text(825,1-(ar/60),area2test{ar},'Color',colorsArea(ar,:),'FontWeight','bold','FontSize',16)
        ylim([chance_l(m)-0.1 1])
        xlim([-150 max(meanPerf{m}(:,2))+100])

    end
    y_ax = [min(keep4legend):(max(keep4legend)-min(keep4legend))/(length(area2test)-1) :max(keep4legend)]
    [aa,bb]=sortrows(keep4legend);
    for ar = 1 : length(area2test)
        text(0,y_ax(ar),area2test{bb(ar)},'Color',colorsArea(bb(ar),:),'FontWeight','bold','FontSize',14,'HorizontalAlignment','right')
    end
    title(['Decoding perf - ' measures{m}])
    xlabel('Nb neurons included')
    ylabel('Average decoding perf')
    set(gca,'FontSize',16)

    subplot(3,6,subplot_idx{m,2});
    for ar = 1 : length(area2test)
        bar(ar,nanmean(param_a{m}(ar,:)),'FaceColor',colorsArea(ar,:));hold on
        line([ar ar],[nanmean(param_a{m}(ar,:))-(nanstd(param_a{m}(ar,:))/sqrt(sum(~isnan(param_a{m}(ar,:))))) nanmean(param_a{m}(ar,:))+(nanstd(param_a{m}(ar,:))/sqrt(sum(~isnan(param_a{m}(ar,:)))))],'Color','k');hold on
    end
    xlim([0 length(area2test)+1])
    set(gca,'Xtick',1:length(area2test),'XtickLabel',area2test,'XtickLabelRotation',30,'FontSize',16)
    title('Info per neuron')

    subplot(3,6,subplot_idx{m,3});
    for ar = 1 : length(area2test)
        bar(ar,nanmean(param_b{m}(ar,:)),'FaceColor',colorsArea(ar,:));hold on
        line([ar ar],[nanmean(param_b{m}(ar,:))-(nanstd(param_b{m}(ar,:))/sqrt(sum(~isnan(param_b{m}(ar,:))))) nanmean(param_b{m}(ar,:))+(nanstd(param_b{m}(ar,:))/sqrt(sum(~isnan(param_b{m}(ar,:)))))],'Color','k');hold on
    end
    xlim([0 length(area2test)+1])
    set(gca,'Xtick',1:length(area2test),'XtickLabel',area2test,'XtickLabelRotation',30,'FontSize',16)
    title('Total info')
end

%% POST HOC - EACH MONKEY - FULL DATASET - STATS

clear

path2go = 'S:\POTT\data\subset-final\'; %- path where SPKpool files are!
stim_1 = load([path2go 'res_LDApop_kfold_final_97sessions.mat'])
stim_2 = load([path2go 'res_LDApop_kfold_final_186sessions.mat'])
stim_all = load([path2go 'res_LDApop_kfold_final_283sessions.mat'])
res_LDApop(1) = stim_1.res_LDApop;
res_LDApop(2) = stim_2.res_LDApop;
res_LDApop(3) = stim_all.res_LDApop;

measures = {'I_chosenjuice' 'I_chosenproba' 'I_chosenside'};
nbNeur = 100;
[colorsArea,colorsArea_sub] = colorMeUp('POTT'); %- COLOR ASSIGNMENT

meanPerf={};
allPerf=cell(3,1);
meanPerf_perm={};
nbNeur_sig =[];
for m = 1 : length(measures)
    x=0;
    for mk = 1 :3
        for ar = 1: length(stim_1.area2test)
            eval(['temp = res_LDApop(mk).' measures{m} '.area_' stim_1.area2test{ar} ';'])
            for n = 1 : length(temp)
                x = x+1;
                pval = sum(temp(n).perf_perm > nanmean(temp(n).perf)) / length(temp(n).perf_perm);
                sem = (nanstd(temp(n).perf)/sqrt(length(temp(n).perf)));
                meanPerf{m}(x,:) = [nanmean(temp(n).perf) temp(n).nbUnit ar nanmean(temp(n).perf)-sem nanmean(temp(n).perf)+sem pval mk];
                sem_perm = (nanstd(temp(n).perf_perm)/sqrt(length(temp(n).perf_perm)));
                meanPerf_perm{m}(x,:) = [nanmean(temp(n).perf_perm) temp(n).nbUnit ar nanmean(temp(n).perf_perm)-sem_perm nanmean(temp(n).perf_perm)+sem_perm mk];
                nbRep = length(temp(n).perf');
                allPerf{m} = [allPerf{m} ; temp(n).perf' repmat(temp(n).nbUnit,nbRep,1) repmat(ar,nbRep,1) repmat(mk,nbRep,1)];
            end
        end
    end
end

for m = 1 : length(measures)
    data_sub = meanPerf{m}(meanPerf{m}(:,end)~=3,[1 2 3 end]);
    mks = {'M' , 'X'}  ;
    %- significance
    modeldata_pop = table(data_sub(:,1),data_sub(:,2),stim_1.area2test(data_sub(:,3))',mks(data_sub(:,4))','VariableNames',{'perf' 'nbunit' 'area' 'mk'});
    models_form = {'perf ~ 1 + area  + (1|mk) + (1|nbunit)' ; 'perf ~ 1 + area  + (1|mk)'};
    % [lme,model_final] = model_comparison(modeldata_pop,models_form,false);
    lme = fitglme(modeldata_pop,models_form{1}); model_final = models_form{1};

    [pval,wald,thr_corr,pval_adj] = area_posthoc(lme,stim_1.area2test,'n');
    disp(['%%%%%%%%%% Decoding perf with model = ' model_final ' %%%%%%%%%%'])
    disp(anova(lme));disp(pval_adj);disp(wald);
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
end

%- Plot the average decoding for each monkey (insets of Fig S1A)
all_decoding_perf = 0;
fig = figure;
for  m = 1  : length(measures)
    subplot(1,5,m);
    decoding_perf = meanPerf{m}(meanPerf{m}(:,2)==nbNeur & meanPerf{m}(:,end)==3,[1 4 5 3]);
    decoding_perf_M = meanPerf{m}(meanPerf{m}(:,2)==nbNeur & meanPerf{m}(:,end)==1,[1 4 5 3]);
    decoding_perf_X = meanPerf{m}(meanPerf{m}(:,2)==nbNeur & meanPerf{m}(:,end)==2,[1 4 5 3]);
    area2plot = stim_1.area2test(decoding_perf(:,end));
    for ar = 1 : length(area2plot)
        if strcmp(area2plot{ar}(1),'1')
            area2plot{ar} = ['a' area2plot{ar}];
        end
    end

    for ar = 1 : length(area2plot)
        bar(decoding_perf(ar,end),decoding_perf(ar,1),'FaceColor',colorsArea(ar,:));hold on
        line([ar ar], [decoding_perf(ar,2) decoding_perf(ar,3) ],'Color','k');hold on
    end
    plot(decoding_perf_M(:,end),decoding_perf_M(:,1),'o','MarkerSize',10,'MarkerFaceColor',[.8 .8 .8],'MarkerEdgeColor','k')
    plot(decoding_perf_X(:,end),decoding_perf_X(:,1),'v','MarkerSize',10,'MarkerFaceColor',[.65 .65 .65],'MarkerEdgeColor','k')
    ylim([0 1])
    xlim([0 length(area2plot)+1])
    set(gca,'Xtick',1:length(area2plot),'XtickLabel',area2plot,'XtickLabelRotation',30,'FontSize',16)
    title(measures{m})
end

%% POSTHOC - SAME BUT FOR SUBSETS NOW

clear
path2go = 'S:\POTT\data\subset-final\'; %- path where SPKpool files are!
load([path2go 'res_LDApopsub_kfold_final_283sessions.mat'])
measures = {'I_chosenside' 'I_unchosenproba' }
%measures = {'P_side' 'P_proba'}
measures = {'I_chosenjuice' 'I_chosenproba' 'I_chosenside'}
nbNeur = 150;

[colorsArea,colorsArea_sub] = colorMeUp('POTT');

clear recall precision specificity f_meas
meanPerf={};
meanPerf_perm={};
nbNeur_sig =[];
for m = 1 : length(measures)
    x=0;
    for ar = 1: length(area2test)
        eval(['temp = res_LDApop.' measures{m} '.area_' area2test{ar} ';'])
        for n = 1 : length(temp)
            if temp(n).nbUnit<=500
            x = x+1;
            pval = sum(temp(n).perf_perm > nanmean(temp(n).perf)) / length(temp(n).perf_perm);
            sem = (nanstd(temp(n).perf)/sqrt(length(temp(n).perf)));

            recall{m}(x,:) = nanmean(cat(1,temp(n).output(:).recall)   );
            precision{m}(x,:) = nanmean(cat(1,temp(n).output(:).precision)   );
            specificity{m}(x,:) = nanmean(cat(1,temp(n).output(:).specificity)   );
            f_meas{m}(x,:) = nanmean(cat(1,temp(n).output(:).f_meas)   );

            meanPerf{m}(x,:) = [nanmean(temp(n).perf) temp(n).nbUnit ar nanmean(temp(n).perf)-sem nanmean(temp(n).perf)+sem pval];
            sem_perm = (nanstd(temp(n).perf_perm)/sqrt(length(temp(n).perf_perm)));
            meanPerf_perm{m}(x,:) = [nanmean(temp(n).perf_perm) temp(n).nbUnit ar nanmean(temp(n).perf_perm)-sem_perm nanmean(temp(n).perf_perm)+sem_perm];
            end
        end
    end
    
    for ar = 1: length(area2test)
        idx = find(meanPerf{m}(:,3)==ar & meanPerf{m}(:,end)<0.05,1,'first');
        if ~isempty(idx)
            nbNeur_sig(m,ar) = meanPerf{m}(idx,2);
        else
            nbNeur_sig(m,ar) = NaN;
        end
        idx2 = find(meanPerf{m}(:,3)==ar ,1,'last');
        nbNeur_test(m,ar) = meanPerf{m}(idx2,2);
    end
end

%- Fig S3C,G
figure;
for m = 1 : 3
    subplot(3,3,[0 3]+m);  chance_l = 0.25;
    for ar = 1: length(area2test)
        dumm = meanPerf{m}(meanPerf{m}(:,3)==ar,[1 2 4 5]);
        plot(dumm(:,2),dumm(:,1),'.-','Color',colorsArea(ar,:),'LineWidth',2,'MarkerSize',15);hold on
        for p = 1:size(dumm,1)
            line([dumm(p,2) dumm(p,2)],[dumm(p,3) dumm(p,4)],'Color',colorsArea(ar,:))
        end
        keep4legend(ar,1) = dumm(1,1);
    
        dumm = meanPerf_perm{m}(meanPerf_perm{m}(:,3)==ar,[1 2 4 5]);
        plot(dumm(:,2),dumm(:,1),'.-','Color',colorsArea(ar,:),'LineWidth',1,'MarkerSize',10);hold on
        for p = 1:size(dumm,1)
            line([dumm(p,2) dumm(p,2)],[dumm(p,3) dumm(p,4)],'Color',colorsArea(ar,:))
        end 
        
        ylim([0.4 1])
        xlim([-150 max(meanPerf{m}(:,2))+100])
    end
    y_ax = [min(keep4legend):(max(keep4legend)-min(keep4legend))/(length(area2test)-1) :max(keep4legend)]
    [aa,bb]=sortrows(keep4legend);
    for ar = 1 : length(area2test)
            text(0,y_ax(ar),area2test{bb(ar)},'Color',colorsArea(bb(ar),:),'FontWeight','bold','FontSize',14,'HorizontalAlignment','right')
    end
    title(['Decoding perf - ' measures{m}])
    xlabel('Nb neurons included')
    ylabel('Average decoding perf')
    set(gca,'FontSize',16)
end

%% STATS FOR SUBSETS AND FIG S3D,H,..

clear

nbNeur = 50;

path2go = 'S:\POTT\data\subset-final\'; %- path where SPKpool files are!
stim_1 = load([path2go 'res_LDApopsub_kfold_final_97sessions.mat'])
stim_2 = load([path2go 'res_LDApopsub_kfold_final_186sessions.mat'])
stim_all = load([path2go 'res_LDApopsub_kfold_final_283sessions.mat'])
res_LDApop(1) = stim_1.res_LDApop;
res_LDApop(2) = stim_2.res_LDApop;
res_LDApop(3) = stim_all.res_LDApop;
measures = {'I_chosenjuice' 'I_chosenproba' 'I_unchosenproba' 'I_chosenside'};

[colorsArea,colorsArea_sub] = colorMeUp('POTT'); %- COLOR ASSIGNMENT

meanPerf={};
allPerf=cell(length(measures),1);
meanPerf_perm={};
nbNeur_sig =[];
for m = 1 : length(measures)
    x=0;
    for mk = 1 :3
        for ar = 1: length(stim_1.area2test)
            eval(['temp = res_LDApop(mk).' measures{m} '.area_' stim_1.area2test{ar} ';'])
            for n = 1 : length(temp)
                x = x+1;
                pval = sum(temp(n).perf_perm > nanmean(temp(n).perf)) / length(temp(n).perf_perm);
                sem = (nanstd(temp(n).perf)/sqrt(length(temp(n).perf)));
                meanPerf{m}(x,:) = [nanmean(temp(n).perf) temp(n).nbUnit ar nanmean(temp(n).perf)-sem nanmean(temp(n).perf)+sem pval mk];
                sem_perm = (nanstd(temp(n).perf_perm)/sqrt(length(temp(n).perf_perm)));
                meanPerf_perm{m}(x,:) = [nanmean(temp(n).perf_perm) temp(n).nbUnit ar nanmean(temp(n).perf_perm)-sem_perm nanmean(temp(n).perf_perm)+sem_perm mk];
                nbRep = length(temp(n).perf');
                allPerf{m} = [allPerf{m} ; temp(n).perf' repmat(temp(n).nbUnit,nbRep,1) repmat(ar,nbRep,1) repmat(mk,nbRep,1)];
            end
        end
    end
end

m=2;
data_sub = meanPerf{m}(meanPerf{m}(:,end)~=3,[1 2 3 end]);
mks = {'M' , 'X'}  ;
%- significance
modeldata_pop = table(data_sub(:,1),data_sub(:,2),stim_1.area2test(data_sub(:,3))',mks(data_sub(:,4))','VariableNames',{'perf' 'nbunit' 'area' 'mk'});
models_form = {'perf ~ 1 + area  + (1|mk) + (1|nbunit)' ; 'perf ~ 1 + area  + (1|mk)'};
% [lme,model_final] = model_comparison(modeldata_pop,models_form,false);
lme = fitglme(modeldata_pop,models_form{1}); model_final = models_form{1};

[pval,wald,thr_corr,pval_adj] = area_posthoc(lme,stim_1.area2test,'n');
disp(['%%%%%%%%%% Decoding perf with model = ' model_final ' %%%%%%%%%%'])
disp(anova(lme));disp(pval_adj);disp(wald);
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

%- Figure S3D,H,L,P
all_decoding_perf = 0;
fig = figure;
for  m = 1  : length(measures)
    subplot(1,5,m);
    decoding_perf = meanPerf{m}(meanPerf{m}(:,2)==nbNeur & meanPerf{m}(:,end)==3,[1 4 5 3]);
    decoding_perf_M = meanPerf{m}(meanPerf{m}(:,2)==nbNeur & meanPerf{m}(:,end)==1,[1 4 5 3]);
    decoding_perf_X = meanPerf{m}(meanPerf{m}(:,2)==nbNeur & meanPerf{m}(:,end)==2,[1 4 5 3]);
    area2plot = stim_1.area2test(decoding_perf(:,end));
    for ar = 1 : length(area2plot)
        if strcmp(area2plot{ar}(1),'1')
            area2plot{ar} = ['a' area2plot{ar}];
        end
    end
    for ar = 1 : length(area2plot)
        bar(decoding_perf(ar,end),decoding_perf(ar,1),'FaceColor',colorsArea(ar,:));hold on
        line([ar ar], [decoding_perf(ar,2) decoding_perf(ar,3) ],'Color','k');hold on
    end
    plot(decoding_perf_M(:,end),decoding_perf_M(:,1),'o','MarkerSize',10,'MarkerFaceColor',[.8 .8 .8],'MarkerEdgeColor','k')
    plot(decoding_perf_X(:,end),decoding_perf_X(:,1),'v','MarkerSize',10,'MarkerFaceColor',[.65 .65 .65],'MarkerEdgeColor','k')
    ylim([0 1])
    xlim([0 length(area2plot)+1])
    set(gca,'Xtick',1:length(area2plot),'XtickLabel',area2plot,'XtickLabelRotation',30,'FontSize',16)
    title(measures{m})
end


