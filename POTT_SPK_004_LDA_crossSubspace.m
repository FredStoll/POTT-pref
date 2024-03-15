%% POTT Pref - Cross subspace decoding of flavor and proba (processing for Figure 4C-D)
%-
%- Author: Fred M. Stoll, Icahn School of Medicine at Mount Sinai, NY
%- Date: 2023.12
%- Related to: Stoll & Rudebeck, Neuron, 2024

clear

path2go = '/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-BehavPrefChange/data-final/neurons/'; %- path where SPKpool files are!
list = dir([path2go '*a_SPKpool.mat']);

param.predic = {'I_juice_by_chosenproba' }; %- select param you want to test..
param.minNeurons = 5; %- min number of neuron to run
param.nComp = 4; %- nb of component to keep for PCA
param.Repetition = 200; %- number of times to run the decoder (100 or 200 ideally)
param.minTrCond = 5; %- minimum number of trials per condition to run
param.overwrite = true;
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

%- check if matrix already exists for each predictor
if exist([path2go 'res_LDA_cross.mat'])==2
    prev = load([path2go 'res_LDA_cross.mat']);
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
            clearvars -except path2go list param areas area2test util s pr res_LDA nSess
            name = list(s).name;

            %- load spike data and histology info
            load([path2go name]); disp(s);
            % load([path2go name(1:end-12) '_AreaInfo.mat']);

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
                if sum(ismember(TrialType_header,param.predic{pr})) ~= 0
                    dumm = TrialType(ismember(TrialType_header,{'nTr' param.predic{pr} }),:)';
                else %- for the proba/juice interaction
                    dumm = TrialType_INS(ismember(TrialType_header,{'nTr' 'I_chosenproba' 'I_chosenjuice' }),:)';
                    dumm = [dumm(:,1) dumm(:,2)+(100*dumm(:,3))];
                end

                %- trials to consider for instrum
                if strcmp(param.predic{pr}(1),'I')
                    when = @(arg,cd) TrialType(strcmp(TrialType_header,arg),:)==cd ;
                    diff_juice = (when('I_juiceL',1) & when('I_juiceR',2)) | (when('I_juiceL',2) & when('I_juiceR',1)); %- take only the different juice trials
                    diff_juice_all = repmat(diff_juice',length(unique(unit_nb)),1);
                end

                predic_all = repmat(dumm(:,2),length(unique(unit_nb)),1);
                factor = [predic_all , unit_nb];

                %- only keep a subset of proba
                remove = [];
                if strcmp(param.predic{pr},'P_proba')
                    keep = [10 30 50 70 90 ];
                    remove = ~ismember(factor(:,1),keep);
                elseif strcmp(param.predic{pr},'I_chosenproba')
                    keep = [30 50 70 90 ];
                    remove = ~ismember(factor(:,1),keep) | ~diff_juice_all;  % was not removing same juice trials before
                elseif strcmp(param.predic{pr},'I_unchosenproba')
                    keep = [10 30 50 70 ];
                    remove = ~ismember(factor(:,1),keep) | ~diff_juice_all;  % was not removing same juice trials before
                elseif strcmp(param.predic{pr},'I_chosenjuice')
                    keep = [1 2];
                    remove = ~diff_juice_all;  % was not removing same juice trials before
                elseif strcmp(param.predic{pr},'I_chosenside')
                    keep = [1 2];
                    remove = ~diff_juice_all;  % was not removing same juice trials before
                elseif strcmp(param.predic{pr},'I_chosenproba_juice') | strcmp(param.predic{pr},'I_chosenproba_by_juice') | strcmp(param.predic{pr},'I_juice_by_chosenproba')
                    keep = [130 150 170 190 230 250 270 290];
                    remove = ~ismember(factor(:,1),keep) | ~diff_juice_all;
                end
                if ~isempty(remove)
                    data(remove,:) = [];
                    factor(remove,:)=[];
                    unit_nb(remove)=[];
                end

                %- mainly used for chosenproba_juice (where both classifiers
                %- should use the same number of trials to avoid biases)
                if ~isempty(factor)
                    tr2cons = factor(factor(:,2)==factor(1,2),1);
                    nCd = length(keep);
                    for k = 1 : nCd
                        nTr(k)=sum(tr2cons==keep(k));
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
                            factor_sub = factor(ismember(unit_nb,takeit),:);

                            if strcmp(param.predic{pr},'I_juice_by_chosenproba') %- in that case, run 2 decoders, 1 per juice!

                                data_centered = NaN(size(data_sub));
                                for nn = 1 : length(takeit)
                                    tempmat = data_sub(factor_sub(:,2)==takeit(nn),:) ;
                                    data_norm = (tempmat -mean(tempmat(:)))/std(tempmat(:));
                                    data_centered(factor_sub(:,2)==takeit(nn),:) = data_norm-mean(data_norm(:));
                                end


                                [XX,YY,param_LDA] = pop_init_noRdmTrials(data_centered,factor_sub,time,'perm',1,'minTr',tr2take,'pop',true); % ,'window',[-.5 1.5]

                                if ~isempty(param.window)
                                    XX_avg{1} = mean(cat(3,XX{time>=param.window(1) & time<=param.window(2)}),3);
                                    XX = XX_avg;
                                    YY = YY(1);
                                end

                                %- proba subspace + decoding flavor
                                clear PCpb
                                allpbs = [30 50 70 90];
                                YY_pb = mod(YY{1},100);
                                for pbs = 1 : length(allpbs)
                                    PCpb.X(pbs,:) = nanmean(XX{1}(YY_pb==allpbs(pbs),:))  ;
                                end
                                PCpb.nComp=1;
                                [PCpb.eigenvectors,PCpb.score,PCpb.eigenvalues,~,PCpb.explained,PCpb.mu] = pca(PCpb.X,'NumComponents',PCpb.nComp);

                                PCpb.Xnew = XX{1};  %- extract new data (all trials and all times)

                                XX_pbproj{1} = (PCpb.Xnew * PCpb.eigenvectors(:,1:PCpb.nComp)) ;
                                YYfl{1} = floor(YY{1}/100);
                                YYpb{1} = mod(YY{1},100);
                                [perf_fl_pbproj(1,p),Post] = lda_kfold(XX_pbproj,YYfl,param);
                                [perf_pb_pbproj(1,p),Post] = lda_kfold(XX_pbproj,YYpb,param);

                                yy_perm = randperm(length(YYfl{1}));
                                for tt = 1 : length(YYfl)
                                    YYperm_fl{tt} = YYfl{1}(yy_perm);
                                end
                                yy_perm = randperm(length(YYpb{1}));
                                for tt = 1 : length(YYpb)
                                    YYperm_pb{tt} = YYpb{1}(yy_perm);
                                end

                                [perf_fl_pbproj_perm(:,p),~] = lda_kfold(XX_pbproj,YYperm_fl,param);
                                [perf_pb_pbproj_perm(:,p),~] = lda_kfold(XX_pbproj,YYperm_pb,param);


                                %- flavor subspace + decoding proba
                                clear PCfl
                                allfl = [1 2];
                                YY_fl = floor(YY{1}/100);
                                for ff = 1 : length(allfl)
                                    PCfl.X(ff,:) = nanmean(XX{1}(YY_fl==allfl(ff),:))  ;
                                end
                                PCfl.nComp=1;
                                [PCfl.eigenvectors,PCfl.score,PCfl.eigenvalues,~,PCfl.explained,PCfl.mu] = pca(PCfl.X,'NumComponents',PCfl.nComp);

                                PCfl.Xnew = XX{1};  %- extract new data (all trials and all times)

                                XX_flproj{1} = (PCfl.Xnew * PCfl.eigenvectors(:,1:PCfl.nComp)) ;
                                YYfl{1} = floor(YY{1}/100);
                                YYpb{1} = mod(YY{1},100);
                                [perf_fl_flproj(1,p),Post] = lda_kfold(XX_flproj,YYfl,param);
                                [perf_pb_flproj(1,p),Post] = lda_kfold(XX_flproj,YYpb,param);


                                [perf_fl_flproj_perm(:,p),~] = lda_kfold(XX_flproj,YYperm_fl,param);
                                [perf_pb_flproj_perm(:,p),~] = lda_kfold(XX_flproj,YYperm_pb,param);

                                angl_PC(p) = rad2deg(subspace(PCfl.eigenvectors(:,1),PCpb.eigenvectors(:,1)));
                            end

                        end

                        if ~strcmp(param.predic{pr},'I_chosenproba_juice')
                            for tt = 1 : size(perf_fl_flproj,1)
                                perf_fl_flproj_pval(tt,1)=sum(perf_fl_flproj_perm(tt,:)>mean(perf_fl_flproj(tt,:)))/size(perf_fl_flproj_perm,2);
                                perf_fl_pbproj_pval(tt,1)=sum(perf_fl_pbproj_perm(tt,:)>mean(perf_fl_pbproj(tt,:)))/size(perf_fl_pbproj_perm,2);
                                perf_pb_flproj_pval(tt,1)=sum(perf_pb_flproj_perm(tt,:)>mean(perf_pb_flproj(tt,:)))/size(perf_pb_flproj_perm,2);
                                perf_pb_pbproj_pval(tt,1)=sum(perf_pb_pbproj_perm(tt,:)>mean(perf_pb_pbproj(tt,:)))/size(perf_pb_pbproj_perm,2);
                            end
                        end
                        %- extract param of interest to save
                        nSess(ar,1) = nSess(ar,1) + 1 ;
                        temp.lda_sess = name(1:8);
                        if ~strcmp(param.predic{pr},'I_chosenproba_juice') & ~strcmp(param.predic{pr},'I_chosenproba_by_juice') & ~strcmp(param.predic{pr},'I_juice_by_chosenproba')
                            temp.Post_all = Post_all;
                            temp.perf = perf;
                            temp.perf_perm = perf_perm;
                            temp.perf_pval = perf_pval;
                        elseif strcmp(param.predic{pr},'I_chosenproba_by_juice') | strcmp(param.predic{pr},'I_juice_by_chosenproba')
                            temp.perf(1,:) = perf_fl_flproj;
                            temp.perf(2,:) = perf_fl_pbproj;
                            temp.perf(3,:) = perf_pb_flproj;
                            temp.perf(4,:) = perf_pb_pbproj;
                            temp.perf_perm = [perf_fl_flproj_perm; perf_fl_pbproj_perm ; perf_pb_flproj_perm;  perf_pb_pbproj_perm];
                            temp.perf_pval = [perf_fl_flproj_pval; perf_fl_pbproj_pval ; perf_pb_flproj_pval;  perf_pb_pbproj_pval];
                        else
                            temp.perf(1,:,:) = perf_j1;
                            temp.perf(2,:,:) = perf_j2;
                        end
                        temp.takeit = takeit;
                        temp.theta = angl_PC;

                        eval(['res_LDA.' param.predic{pr} '.' area2test{ar} '(nSess(ar)) = temp;'])
                    end
                end
            end
        end
        %- add param
        eval(['res_LDA.' param.predic{pr} '.param = param;'])
    end
    %- final save!
    save([path2go 'res_LDA_cross.mat'],'res_LDA','list','area2test','path2go')
end

