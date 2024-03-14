%% POTT Pref - Population decoding and preference modulation - Figure 8
%-
%- Require the time windows for PRE/POST pref change extracted in
%- POTT_BHV_001_GLM (Matrix Pref_bins.mat)
%-
%- Author: Fred M. Stoll, Icahn School of Medicine at Mount Sinai, NY
%- Date: 2023.08
%- Related to: Stoll & Rudebeck, Neuron, 2024

clear

path2go = '/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-BehavPrefChange/data/neurons/subset-final/'; %- path where SPKpool files are!
list = dir([path2go '*a_SPKpool.mat']);

param.predic = { 'I_chosenproba_pref' 'I_chosenjuice_pref' }; %- select param you want to test..
param.minNeurons = 5; %- min number of neuron to run
param.nComp = 4; %- nb of component to keep for PCA
param.Repetition = 200; %- number of times to run the decoder (100 or 200 ideally)
param.minTrCond = 5; %- minimum number of trials per condition to run
param.overwrite = false;
param.thresh = false;
param.thr = 1; % in Hz, used if param.thresh is 'y'
param.thr_len = 500; % in ms, used if param.thresh is 'y'
param.rmv = 100; % remove xxx ms on each side on every events (avoid overlaps between bins and smoothing problems)
param.bins4decoding = [2] ; %- perform decoding on subset on bins (around stim period here)
param.window = [200 800]; %- empty bracket if every time points, otherwise will do LDA on average FR during [t_start t_end] in ms
param.normalize_fr = false;

disp(['Computing LDA on ' num2str(length(list)) ' sessions'])
areas = utils_POTT_areas;

area2test = {'vlPFC' 'OFC' 'IFG' 'LAI' 'AMG'};

util.minmaxnorm = @(data) (data-min(data))/(max(data)-min(data));
util.meannorm  = @(data) (data-mean(data))/(max(data)-min(data));
util.minnorm  = @(data) (data-min(min(data)))/(max(max(data))-min(min(data)));

%- load preference periods
path2go2 = 'C:\\Users\fred\Dropbox\Rudebeck Lab\ANA-POTT-BehavPrefChange\data\';
load([path2go2 'Pref_bins.mat'])
pref_ch_mk = [pref_ch_mk{1} ;pref_ch_mk{2}];

%- check if matrix already exists for each predictor
if exist([path2go 'res_LDA_stim_pref.mat'])==2
    prev = load([path2go 'res_LDA_stim_pref.mat']);
    done = [];
    for pr = 1 : length(param.predic)
        if isfield(prev.res_LDA,param.predic{pr})
            done(pr) = true;
        end
    end
    if ~param.overwrite
        param.predic(find(done==true))=[];
    end
    res_LDA = prev.res_LDA;
end
if param.overwrite
    clear res_LDA
end

if ~isempty(param.predic) %- if there is some predictors not done yet

    for pr = 1 : length(param.predic)
        nSess = zeros(length(area2test),1);
        for s = 1 : length(list)
            clearvars -except path2go list param areas area2test util s pr res_LDA nSess pref_ch_mk
            name = list(s).name;

            if sum(ismember(pref_ch_mk(:,1),name(1:8)))==1
                sessnum = find(ismember(pref_ch_mk(:,1),name(1:8)));
                pref_curr = pref_ch_mk(sessnum,:);

                %- load spike data and histology info
                load([path2go name]); disp(s);

                if ~isempty(neurons_area) %- skip sessions where no recorded neurons where in area2test!

                    %- create time vector and remove overlapping windows/border events
                    bins2remove = (param.rmv/subsp)-1;

                    n_evt = length(times_evts);
                    lin=length(neurons_rank)*n_evt;
                    nTimes = [0 ; (sum(abs(times(:,:)),2)/subsp)];
                    col = sum(nTimes);
                    bins = NaN(1,col);
                    time = NaN(1,col);

                    for b = 1 : length(times_evts)
                        time(sum(nTimes(1:b))+1:sum(nTimes(1:b+1))) =  (times(b,1):subsp:times(b,2)-(subsp/2));
                        bins(sum(nTimes(1:b))+1:sum(nTimes(1:b+1))) = b;
                        bins([sum(nTimes(1:b))+1:sum(nTimes(1:b))+1+bins2remove   sum(nTimes(1:b+1))-bins2remove:sum(nTimes(1:b+1))    ]) = NaN; %- remove 100 ms each side (avoid overlaps between bins and smoothing problems)
                    end

                    %- extract the data for the task 2 consider
                    if strcmp(param.predic{pr}(1),'I') %- instrumental
                        SPK_data = SPK_INS;
                        Tr_Clust_data = Tr_Clust_INS;
                        TrialType = TrialType_INS;
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
                        for n = 1 : size(neurons_rank,1)
                            if isempty(findenough(mean(SPK_data_cut(Tr_Clust_data(:,2)==n,:)),param.thr,param.thr_len/subsp,'>='))
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
                    unit_nb = Tr_Clust_data(:,2);

                    dumm = TrialType_INS(ismember(TrialType_header,{'nTr' 'I_chosenproba' 'I_chosenjuice' }),:)';
                    dumm = [dumm(:,1) dumm(:,2)+(100*dumm(:,3))];

                    when = @(arg,cd) TrialType_INS(strcmp(TrialType_header,arg),:)==cd ;
                    diff_juice = (when('I_juiceL',1) & when('I_juiceR',2)) | (when('I_juiceL',2) & when('I_juiceR',1)); %- take only the different juice trials

                    %- extract only the diff juice trials
                    tt = find(diff_juice);

                    keep = [130 150 170 190 230 250 270 290 ];
                    %  remove = ~ismember(dumm(:,2),keep) | ~diff_juice';

                    preTr = pref_curr{1,2}; %- out of the diff juice trials
                    preTr4spk = tt(preTr); %- trial to take for neurons PRE
                    preTr_fact = dumm(preTr4spk,:);
                    preRemove = ~ismember(preTr_fact(:,2),keep) ;
                    preTr_fact(preRemove,:)=[];

                    postTr = pref_curr{1,3}; %- out of the diff juice trials
                    postTr4spk = tt(postTr); %- trial to take for neurons POST
                    postTr_fact = dumm(postTr4spk,:);
                    postRemove = ~ismember(postTr_fact(:,2),keep) ;
                    postTr_fact(postRemove,:)=[];

                    if pref_curr{1,4}(2)-pref_curr{1,4}(1)>0
                        allfactors = [preTr_fact(:,1) mod(preTr_fact(:,2),100) abs(floor(preTr_fact(:,2)/100)-3) ones(size(preTr_fact(:,1))) ; ...
                            postTr_fact(:,1) mod(postTr_fact(:,2),100) abs(floor(postTr_fact(:,2)/100)-3) 2*ones(size(postTr_fact(:,1)))]; %- inverse juice (1 = valued / 2 = devalued)
                    else
                        allfactors = [preTr_fact(:,1) mod(preTr_fact(:,2),100) floor(preTr_fact(:,2)/100) ones(size(preTr_fact(:,1))) ; ...
                            postTr_fact(:,1) mod(postTr_fact(:,2),100) floor(postTr_fact(:,2)/100) 2*ones(size(postTr_fact(:,1)))]; %- no inverse juice (1 = valued / 2 = devalued)
                    end

                    data = data( ismember(Tr_Clust_INS(:,1),allfactors(:,1)),:);
                    unit_nb = unit_nb(ismember(Tr_Clust_INS(:,1),allfactors(:,1)));
                    factor = repmat(allfactors,length(neurons_area),1);

                    if strcmp(param.predic{pr},'I_chosenjuice_pref')
                        factor_comb = [1000*factor(:,4) + 100*factor(:,3) unit_nb];
                        all_keep = [1100 1200 2100 2200];
                    elseif strcmp(param.predic{pr},'I_chosenproba_pref')
                        factor_comb = [1000*factor(:,4) + factor(:,2)  unit_nb];
                        all_keep = [1030 1050 1070 1090 2030 2050 2070 2090];
                    end

                    %- mainly used for chosenproba_juice (where both classifiers
                    %- should use the same number of trials to avoid biases)
                    if ~isempty(factor)
                        tr2cons = factor_comb(unit_nb==unit_nb(1,1),1);
                        nCd = length(all_keep);nTr = [];
                        for k = 1 : nCd
                            nTr(k)=sum(tr2cons==all_keep(k));
                        end
                        tr2take = min(nTr);
                    end
                    all_units = unique(unit_nb);

                    %- for every area with enough neurons... min is 5
                    for ar = 1 : length(area2test)
                        eval(['takeit = all_units(find(ismember(neurons_area,areas.' area2test{ar} ')));'])

                        if length(takeit) >= param.minNeurons & tr2take >= param.minTrCond %- if enough neurons + trials

                            %- Initialize
                            if ~isempty(param.window)
                                Post_all = zeros(nCd,nCd);
                            else
                                Post_all = zeros(size(data,2),nCd,nCd);
                            end
                            for p = 1 : param.Repetition
                                disp(['session ' num2str(s) '/' num2str(length(list)) ' - ' area2test{ar} ' - ' param.predic{pr} ' - perm ' num2str(p) '/' num2str(param.Repetition)])

                                %- take neurons from the area 2 test
                                data_sub = data(ismember(unit_nb,takeit),:);
                                factor_sub = factor_comb(ismember(unit_nb,takeit),:);

                                if ~strcmp(param.predic{pr},'I_chosenproba_juice') & ~strcmp(param.predic{pr},'I_chosenproba_pref') & ~strcmp(param.predic{pr},'I_chosenjuice_pref')
                                    %- reformat + use only a similar number of trials for each condition
                                    [XX,YY,param_LDA] = pop_init_noRdmTrials(data_sub,factor_sub,time,'perm',1,'minTr',tr2take,'pop',true); % ,'window',[-.5 1.5]

                                    if ~isempty(param.window)
                                        XX_avg{1} = mean(cat(3,XX{time>=param.window(1) & time<=param.window(2)}),3);
                                        XX = XX_avg;
                                        YY = YY(1);
                                    end
                                    [perf(:,p),Post] = pca_lda_kfold(XX,YY,param);

                                    Post_all =  Post_all + Post;


                                    yy_perm = randperm(length(YY{1}));
                                    for tt = 1 : length(YY)
                                        YYperm{tt} = YY{1}(yy_perm);
                                    end
                                    [perf_perm(:,p),~] = pca_lda_kfold(XX,YYperm,param);

                                else %- in that case, run 2 decoders, 1 per juice!

                                    %- decoder juice 1
                                    data_sub_j1 = data_sub;
                                    factor_sub_j1 = factor_sub;
                                    if strcmp(param.predic{pr},'I_chosenjuice_pref')
                                        j1 = [1100 1200];
                                    elseif strcmp(param.predic{pr},'I_chosenproba_pref')
                                        j1 = [1030 1050 1070 1090];
                                    end
                                    remove = ~ismember(factor_sub_j1(:,1),j1);
                                    data_sub_j1(remove,:) = [];
                                    factor_sub_j1(remove,:)=[];

                                    [XX,YY,param_LDA] = pop_init_noRdmTrials(data_sub_j1,factor_sub_j1,time,'perm',1,'minTr',tr2take,'pop',true); % ,'window',[-.5 1.5]

                                    if ~isempty(param.window)
                                        XX_avg{1} = mean(cat(3,XX{time>=param.window(1) & time<=param.window(2)}),3);
                                        XX = XX_avg;
                                        YY = YY(1);
                                    end
                                    [perf_j1(:,p),~] = pca_lda_kfold(XX,YY,param);

                                    %- decoder juice 2
                                    data_sub_j2 = data_sub;
                                    factor_sub_j2 = factor_sub;
                                    if strcmp(param.predic{pr},'I_chosenjuice_pref')
                                        j2 = [2100 2200];
                                    elseif strcmp(param.predic{pr},'I_chosenproba_pref')
                                        j2 = [2030 2050 2070 2090];
                                    end
                                    remove = ~ismember(factor_sub_j2(:,1),j2);
                                    data_sub_j2(remove,:) = [];
                                    factor_sub_j2(remove,:)=[];

                                    [XX,YY,param_LDA] = pop_init_noRdmTrials(data_sub_j2,factor_sub_j2,time,'perm',1,'minTr',tr2take,'pop',true); % ,'window',[-.5 1.5]

                                    if ~isempty(param.window)
                                        XX_avg{1} = mean(cat(3,XX{time>=param.window(1) & time<=param.window(2)}),3);
                                        XX = XX_avg;
                                        YY = YY(1);
                                    end
                                    [perf_j2(:,p),~] = pca_lda_kfold(XX,YY,param);
                                end


                            end

                            if ~strcmp(param.predic{pr},'I_chosenproba_juice') & ~strcmp(param.predic{pr},'I_chosenproba_pref') & ~strcmp(param.predic{pr},'I_chosenjuice_pref')
                                for tt = 1 : size(perf,1)
                                    perf_pval(tt,1)=sum(perf_perm(tt,:)>mean(perf(tt,:)))/size(perf_perm,2);
                                end
                            end

                            %- extract param of interest to save
                            nSess(ar,1) = nSess(ar,1) + 1 ;
                            temp.lda_sess = name(1:8);
                            if ~strcmp(param.predic{pr},'I_chosenproba_juice') & ~strcmp(param.predic{pr},'I_chosenproba_pref') & ~strcmp(param.predic{pr},'I_chosenjuice_pref')
                                temp.Post_all = Post_all;
                                temp.perf = perf;
                                temp.perf_perm = perf_perm;
                                temp.perf_pval = perf_pval;
                            else
                                temp.perf(1,:,:) = perf_j1;
                                temp.perf(2,:,:) = perf_j2;
                                temp.pref_dir = pref_curr{4};
                                temp.rt_dir = pref_curr{5};
                                temp.ft_dir = pref_curr{6};
                            end
                            temp.takeit = takeit;
                            eval(['res_LDA.' param.predic{pr} '.' area2test{ar} '(nSess(ar)) = temp;'])

                        end
                    end
                end

            end
        end
        %- add param
        eval(['res_LDA.' param.predic{pr} '.param = param;'])
    end
    %- final save!
    save([path2go 'res_LDA_stim_pref.mat'],'res_LDA','list','area2test','path2go')
end


%% Post hoc!

clear

path2go = '/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-BehavPrefChange/data/neurons/subset-final/'

%- load the decoding results
measures = {'I_chosenproba_pref'  'I_chosenjuice_pref' };
name = {'perf_pbxju' 'perf_pb' }
load([path2go 'res_LDA_stim_pref.mat']);

area2test = {'vlPFC' 'OFC' 'IFG' 'LAI' 'AMG'}
[colorsArea,colorsArea_sub] = colorMeUp('POTT'); %- COLOR ASSIGNMENT

maxPerf = [];
sessions_all=char();
for m = 1 : length(measures)
    tab(m).modeldata=[];

    for ar = 1 : length(area2test)
        eval(['curr = res_LDA.' measures{m} '.' area2test{ar} ';'])
        ss = size(curr(1).perf);
        perf = [];
        pref=[];
        sess=[];
        mk = [];
        perf_raw_pval=[];
        for s = 1 : length(curr)
            perf(s,:) = nanmean(curr(s).perf,length(ss))';
            pref(s,:) = curr(s).pref_dir;

            sess{s,1} = curr(s).lda_sess ;
            mk{s,1} = curr(s).lda_sess(1) ;
        end
            curr_area = repmat(area2test(ar),size(mk));

        tab(m).modeldata = [tab(m).modeldata ; table(perf(:,1),perf(:,2),pref(:,1),pref(:,2),curr_area,sess,mk,'VariableNames',{'perf1' 'perf2' 'pref1' 'pref2' 'area' 'session' 'mk'})];
    end
end

%- Check relation between proba and juice (Fig 8D)
scatt_pbju=[];scatt_area={};
x=0;
for i = 1 : height(tab(2).modeldata)
    findit = ismember(tab(1).modeldata.session,tab(2).modeldata.session{i}) & ismember(tab(1).modeldata.area,tab(2).modeldata.area{i});
    if sum(findit==1)
        x=x+1;
        if find(abs([tab(2).modeldata.pref1(i) tab(2).modeldata.pref2(i)]-0.5)==max(abs([tab(2).modeldata.pref1(i) tab(2).modeldata.pref2(i)]-0.5)))==2 %- increase in pref
            scatt_pbju(x,:)=[tab(1).modeldata.perf2(findit)-tab(1).modeldata.perf1(findit) tab(2).modeldata.perf2(i)-tab(2).modeldata.perf1(i)];
        else
          scatt_pbju(x,:)=[tab(1).modeldata.perf1(findit)-tab(1).modeldata.perf2(findit) tab(2).modeldata.perf1(i)-tab(2).modeldata.perf2(i)];
        end
        scatt_area{x}= tab(2).modeldata.area{i};
    end
end

figure
for ar = 1 : length(area2test)-1
    subplot(2,2,ar)

    scatt_pbju_sub = scatt_pbju(ismember(scatt_area,area2test(ar)),:);
    scatt_pbju_sub(sum(isnan(scatt_pbju_sub),2)~=0,:)=[];
    plot(scatt_pbju_sub(:,1),scatt_pbju_sub(:,2),'.','Color',colorsArea(ar,:));hold on
    [r,p]=corrcoef(scatt_pbju_sub(:,1),scatt_pbju_sub(:,2));
    length(scatt_pbju_sub(:,1))
     cf = fit(scatt_pbju_sub(:,1),scatt_pbju_sub(:,2),'poly1'); % fit
     p_handle = plot(cf,'k','predfunc',.95); % Plot fit
        set(p_handle,'Color',colorsArea(ar,:),'LineWidth',2);
        legend off
    title([area2test{ar} ' r=' num2str(r(2,1)) ,' p=' num2str(p(2,1))])
    line([-0.3 0.3],[0 0],'Color','k');hold on
        line([0 0],[-0.3 0.3],'Color','k')
    xlim([-0.3 0.3])
    ylim([-0.3 0.3]);box on
    xlabel('ChosenProba difference')
    ylabel('ChosenJuice difference')
end

%- categories of effects (increase proba / decrease flavor ...)
nCat=[];
for ar = 1 : length(area2test)-1
    scatt_pbju_sub = scatt_pbju(ismember(scatt_area,area2test(ar)),:);
    scatt_pbju_sub(sum(isnan(scatt_pbju_sub),2)~=0,:)=[];

   nCat(ar,:) = [ sum(scatt_pbju_sub(:,1)<0 & scatt_pbju_sub(:,2)<0) , ...
    sum(scatt_pbju_sub(:,1)<0 & scatt_pbju_sub(:,2)>0) ,...
    sum(scatt_pbju_sub(:,1)>0 & scatt_pbju_sub(:,2)>0) ,...
    sum(scatt_pbju_sub(:,1)>0 & scatt_pbju_sub(:,2)<0) , ...
    length(scatt_pbju_sub(:,1))];
end

pval_area =[]; chi2stat_area=[];
for ar = 1 : length(area2test)-1
    [~,chi2stat_area(ar,1),pval_area(ar,1)] = chi2_fms(nCat(ar,1),nCat(ar,end),nCat(ar,2),nCat(ar,end));
    [~,chi2stat_area(ar,2),pval_area(ar,2)] = chi2_fms(nCat(ar,1),nCat(ar,end),nCat(ar,3),nCat(ar,end));
    [~,chi2stat_area(ar,3),pval_area(ar,3)] = chi2_fms(nCat(ar,1),nCat(ar,end),nCat(ar,4),nCat(ar,end));
    [~,chi2stat_area(ar,4),pval_area(ar,4)] = chi2_fms(nCat(ar,2),nCat(ar,end),nCat(ar,3),nCat(ar,end));
    [~,chi2stat_area(ar,5),pval_area(ar,5)] = chi2_fms(nCat(ar,2),nCat(ar,end),nCat(ar,4),nCat(ar,end));
    [~,chi2stat_area(ar,6),pval_area(ar,6)] = chi2_fms(nCat(ar,3),nCat(ar,end),nCat(ar,4),nCat(ar,end));
    [h, crit_p, pval_area(ar,:)]=fdr_bh(pval_area(ar,:),.05);

end

figure;
for ar = 1:length(area2test)-1
    subplot(4,1,ar)
    b = bar(nCat(ar,1:end-1)./nCat(ar,end));
    b.FaceColor = colorsArea(ar,:);
    ylim([0 0.5])
end

%- count number for table
for ar = 1 :length(area2test)
    nbSess(ar,1) = length(unique(tab(2).modeldata.session(ismember(tab(2).modeldata.area,area2test(ar)) &  ismember(tab(2).modeldata.mk ,'M')       )));
    nbSess(ar,2) = length(unique(tab(2).modeldata.session(ismember(tab(2).modeldata.area,area2test(ar)) &  ismember(tab(2).modeldata.mk ,'X')       )));
end

%- angle plot (Figure 8A insets)
figure;
x=0;pval=[];z=[];
for m = 1 : length(measures)
    for ar = 1 : length(area2test)-1
        x=x+1;
        modeldata = tab(m).modeldata(ismember(tab(m).modeldata.area,area2test{ar}),:);
        perfCh=[];prefCh=[];
        for i = 1 : height(modeldata)
            if find(abs([modeldata.pref1(i) modeldata.pref2(i)]-0.5)==max(abs([modeldata.pref1(i) modeldata.pref2(i)]-0.5)))==2 %- increase in pref
               perfCh(i) = modeldata.perf2(i)-modeldata.perf1(i);
               prefCh(i) = abs(modeldata.pref2(i)-.5)-abs(modeldata.pref1(i)-.5);
            else
               perfCh(i) =  modeldata.perf1(i)-modeldata.perf2(i);
               prefCh(i) = abs(modeldata.pref1(i)-.5)-abs(modeldata.pref2(i)-.5);
            end
          
        end
        subplot(2,length(area2test)-1,x)
        circ_plot(atan2(perfCh, prefCh),'hist',[],20,true,true,'linewidth',2,'color',colorsArea(ar,:))
        [pval(m,ar)]= circ_symtest(atan2(perfCh, prefCh));
    end
end

%- Figure 8BC - could make it prettier.........!
chi2stat=[];pval=[];pCh=[];
exp_dir = {'left' 'right'};
for m = 1 : length(measures)
    figure;
    for ar = 1 : length(area2test)-1
        modeldata = tab(m).modeldata(ismember(tab(m).modeldata.area,area2test{ar}),:);
        prefCh=[];
        perfCh=[];
        perfCh_both=[];
        for i = 1 : height(modeldata)
            if find(abs([modeldata.pref1(i) modeldata.pref2(i)]-0.5)==max(abs([modeldata.pref1(i) modeldata.pref2(i)]-0.5)))==2 %- increase in pref
               perfCh(i) = modeldata.perf2(i)-modeldata.perf1(i);
               perfCh_both(i,:) = [modeldata.perf2(i) modeldata.perf1(i)];
               prefCh(i) = abs(modeldata.pref2(i)-.5)-abs(modeldata.pref1(i)-.5);
            else
               perfCh(i) =  modeldata.perf1(i)-modeldata.perf2(i);
               perfCh_both(i,:) = [modeldata.perf1(i) modeldata.perf2(i)];
               prefCh(i) = abs(modeldata.pref1(i)-.5)-abs(modeldata.pref2(i)-.5);
            end
        end
        pCh(ar,:)=[sum(sign(perfCh)==-1) sum(sign(perfCh)==1)];
         [f,ep]= ksdensity(perfCh);
         plot(ep,f,'Color',colorsArea(ar,:));
        %plot(prefCh,perfCh,'o','Color',colorsArea(ar,:));
        hold on
        [tbl,chi2stat(ar,m),pval(ar,m)] = chi2_fms(pCh(ar,1),sum(pCh(ar,:)),pCh(ar,2),sum(pCh(ar,:)));
        [p_Wilcox(ar,m),h,sta_Wilcox(ar,m)]=signrank(perfCh_both(:,1),perfCh_both(:,2));
        Z_Wilcox(ar,m) = sta_Wilcox(ar,m).zval;
    end
    xline(0)
    disp(pCh./repmat(sum(pCh')',1,2))
    figure;
    barh(pCh./repmat(sum(pCh')',1,2),'stacked')
    axis ij
end

%- Figure 8A
colorsPosNeg = cbrewer('qual', 'Set1', 3);
figure;
x = 0;
for m = 1 : length(measures)
    for ar = 1 : length(area2test)-1
        modeldata = tab(m).modeldata(ismember(tab(m).modeldata.area,area2test{ar}),:);
        x = x + 1 ;
        subplot(2,length(area2test)-1,x)

        for i = 1 : height(modeldata)
            p1 = [abs(modeldata.pref1(i)-.5) modeldata.perf1(i)];
            p2 = [abs(modeldata.pref2(i)-.5) modeldata.perf2(i)];
            dp = p2-p1;    
            if sum(sign(dp))==0
                plot([p1(1) p2(1)],[p1(2) p2(2)],'Color',colorsPosNeg(1,:));hold on
                plot(p2(1),p2(2),'.','MarkerSize',10,'Color',colorsPosNeg(1,:))
                % quiver(p1(1),p1(2),dp(1),dp(2),0,'Color',colorsPosNeg(1,:));hold on
            else
                plot([p1(1) p2(1)],[p1(2) p2(2)],'Color',colorsPosNeg(2,:));hold on
                plot(p2(1),p2(2),'.','MarkerSize',10,'Color',colorsPosNeg(2,:))
                % quiver(p1(1),p1(2),dp(1),dp(2),0,'Color',colorsPosNeg(2,:));hold on
            end
        end
        if m == 1 ; ylim([.15 .9])
        else ; ylim([.4 1])
        end
        xlim([0 .5])
        title(area2test{ar})
    end
end
