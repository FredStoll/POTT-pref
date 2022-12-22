
%% Pref ANOVA - neurons


clear

monkey = 'X'
norm_me = @(data) -1+((data-min(data))*2)/(max(data)-min(data)) ;
norm_pb = @(data) -1+((data-0.1)*2)/(0.9-0.1) ;
predic_INS={'I_chosenproba'; 'I_unchosenproba' ;'I_chosenjuice'}; %- select any param you want to test.. order matters!

% path2go = '/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-BehavPrefChange/data/';
% %path2go = 'C:\\Users\fred\Dropbox\Rudebeck Lab\ANA-POTT-BehavPrefChange\data\';
% load([path2go 'Mimic_behav_norm_ALL_prevJuice_only.mat'])

path2go = '/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-BehavPrefChange/data/neurons/subset/'; %- path where SPKpool files are!
path2go2 = '/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-BehavPrefChange/data/';
load([path2go2 'Pref_bins.mat'])

if strcmp(monkey,'X')
    load([path2go2 'Mimic_behav_bins.mat'])
 %   load([path2go2 'Mimic_behav_bins_perms.mat'],'ALL','param')
    list = dir([path2go 'X*a_SPKpool.mat']);
    pref_ch = pref_ch_mk{2};
else
    load([path2go2 'Morbier_behav_bins.mat'])
 %   load([path2go2 'Morbier_behav_bins_perms.mat'],'ALL','param')
    list = dir([path2go 'M*a_SPKpool.mat']);
    pref_ch = pref_ch_mk{1};
end

clear days
for i = 1 : length(ALL)
    days{i,:} = [ALL(i).name(end-17:end-11) 'a'];
end

x = 0;

rmv = 100; % remove xxx ms on each side on every events (avoid overlaps between bins and smoothing problems)
% times_evts = {'FixFP_onset' 'Stim_onset' 'Resp_onset' 'FixResp_onset' 'FB_onset' 'Rew_onset' 'FB_offset'};
bins4decoding=[2]; %- perform decoding on subset on bins (stim and go period here)
binwindow = [200 800];
area_list = utils_POTT_areas;
area2test = {'a11ml' 'a13ml' 'a12mr' 'a12ol' 'IFG' 'LAI' 'AMG' 'Cd' 'PUT' }; % change order
area2test = {'vlPFC' 'OFC' 'IFG' 'LAI' 'Striatum' 'AMG'};

minNeurons = 5;

all_unit = cell(length(area2test),length(list),3);
all_unit_ft = cell(length(area2test),length(list),3);
all_unit_rt = cell(length(area2test),length(list),3);
pref_sig = NaN(length(list),4);
pca_explained = cell(length(area2test),length(list),3);
all_unit_data=cell(length(area2test),length(list),3);

meannorm  = @(data) (data-mean(mean(data)))/(max(max(data))-min(min(data)));
minnorm  = @(data) (data-min(min(data)))/(max(max(data))-min(min(data)));
normal = true;

final_data=table();

measures = {'all'} ;
clear res_anova 
ne = 0;
for s =  1 : length(list)
    disp(s)
    SPK = load([path2go list(s).name]);
    idxSess = ismember(days,list(s).name(1:8));
    takebhv = ismember(pref_ch(:,1),list(s).name(1:8));
    if ~isempty(SPK.neurons_area) && sum(takebhv)==1  %- skip sessions where no recorded neurons where in area2test!
        % spkfile = [ ALL(s).name(end-17:end-10) '_SPKpool.mat']
        % areafile = [ ALL(s).name(end-17:end-10) '_AreaInfo.mat']
        
        %load([path2go '/neurons/M072718a_SPKpool.mat'])
        
        %- find times
        bins2remove = (rmv/SPK.subsp)-1;
        
        n_evt = length(SPK.times_evts);
        lin=length(SPK.neurons_rank)*n_evt;
        nTimes = [0 ; (sum(abs(SPK.times(:,:)),2)/SPK.subsp)];
        col = sum(nTimes);
        bins = NaN(1,col);
        time = NaN(1,col);
        
        for b = 1 : length(SPK.times_evts)
            time(sum(nTimes(1:b))+1:sum(nTimes(1:b+1))) =  (SPK.times(b,1):SPK.subsp:SPK.times(b,2)-(SPK.subsp/2));
            bins(sum(nTimes(1:b))+1:sum(nTimes(1:b+1))) = b;
            bins([sum(nTimes(1:b))+1:sum(nTimes(1:b))+1+bins2remove   sum(nTimes(1:b+1))-bins2remove:sum(nTimes(1:b+1))    ]) = NaN; %- remove 100 ms each side (avoid overlaps between bins and smoothing problems)
        end
        
        %- keep only considered bins
        % SPK.SPK_INS(:,~ismember(bins,bins4decoding))=[];
        time2cons = ismember(bins,bins4decoding) & time>=binwindow(1) & time<=binwindow(2);
        SPK.SPK_INS = SPK.SPK_INS(:,time2cons);
        time(~time2cons) = [];
        
        
        %  SPK.SPK_INS = SPK.SPK_INS(:,296:316);
        
        
        %- find trial ID included for each 30-trial behavioral model
%         clear trID
%         dumm = find(ALL(idxSess).trial2take);
%         for b = 1 : length(dumm)-param.binSize
%             trID(b,:) = dumm(b:b+param.binSize);
%         end
        
        % find trial ID for the neurons
        trID_neurons = find(ALL(idxSess).TrialType(1,:)==2 & ALL(idxSess).TrialType(2,:)==0);
        
        %     figure;plot(ALL(s).TrialType(11,trID_neurons));hold on
        %     plot(SPK.TrialType_INS(12,:))
        [ALL(idxSess).Z_trend_ft,~]=Mann_Kendall(ALL(idxSess).ft_bins,0.01);
        [ALL(idxSess).Z_trend_rt,~]=Mann_Kendall(ALL(idxSess).rt_bins,0.01);
        [ALL(idxSess).Z_trend_pref,~]=Mann_Kendall(ALL(idxSess).pref_bins,0.01);
        
        if length(trID_neurons)~=size(SPK.TrialType_INS,2)
            warning('Problem in # of trials......')
            x = x + 1;
            pb_files(x)=s;
            
        else
            
            pref_ch(takebhv,:)
            
            dumm = SPK.TrialType_INS(ismember(SPK.TrialType_header,{'nTr' 'I_chosenproba' 'I_chosenjuice' }),:)';
            dumm = [dumm(:,1) dumm(:,2)+(100*dumm(:,3))];
            
            when = @(arg,cd) SPK.TrialType_INS(strcmp(SPK.TrialType_header,arg),:)==cd ;
            diff_juice = (when('I_juiceL',1) & when('I_juiceR',2)) | (when('I_juiceL',2) & when('I_juiceR',1)); %- take only the different juice trials
            
            %- extract only the diff juice trials
            tt = find(diff_juice);
            
            keep = [130 150 170 190 230 250 270 290];
          %  remove = ~ismember(dumm(:,2),keep) | ~diff_juice';
            
            preTr = pref_ch{takebhv,2}; %- out of the diff juice trials
            preTr4spk = tt(preTr); %- trial to take for neurons
            preTr_fact = dumm(preTr4spk,:);
            preRemove = ~ismember(preTr_fact(:,2),keep) ;
            preTr_fact(preRemove,:)=[];
          
            postTr = pref_ch{takebhv,3}; %- out of the diff juice trials
            postTr4spk = tt(postTr); %- trial to take for neurons
            postTr_fact = dumm(postTr4spk,:);
            postRemove = ~ismember(postTr_fact(:,2),keep) ;
            postTr_fact(postRemove,:)=[];
            
            if pref_ch{takebhv,4}(2)-pref_ch{takebhv,4}(1)>0
                allfactors = [preTr_fact(:,1) norm_pb(mod(preTr_fact(:,2),100)/100) abs(floor(preTr_fact(:,2)/100)-3) ones(size(preTr_fact(:,1))) ; ...
                              postTr_fact(:,1) norm_pb(mod(postTr_fact(:,2),100)/100) abs(floor(postTr_fact(:,2)/100)-3) 2*ones(size(postTr_fact(:,1)))]; %- inverse juice (1 = valued / 2 = devalued)
            else
                allfactors = [preTr_fact(:,1) norm_pb(mod(preTr_fact(:,2),100)/100) floor(preTr_fact(:,2)/100) ones(size(preTr_fact(:,1))) ; ...
                              postTr_fact(:,1) norm_pb(mod(postTr_fact(:,2),100)/100) floor(postTr_fact(:,2)/100) 2*ones(size(postTr_fact(:,1)))]; %- no inverse juice (1 = valued / 2 = devalued)
            end
            

            clear Ap PEVs Fs Omega2 R2 Betas
            for n = 1 : length(SPK.neurons_area)
                
                data = mean(SPK.SPK_INS(SPK.Tr_Clust_INS(:,2)==n,:),2);
                
                [Ap(:,n),At,Amdl,~] = anovan(data(allfactors(:,1)),allfactors(:,2:4),'continuous',[1],'varname',{'ChosenProba' 'Val_Deval' 'Pre_Post'},'model',[1 0 0;0 1 0;0 0 1; 1 1 0 ; 1 0 1 ; 0 1 1 ; 1 1 1],'display','off');
                for m = 1 : length(Ap(:,n))
                    PEVs(m,n) = ((At{m+1,2}/At{end,2}))*100;
                    Fs(m,n) = At{m+1,6};
                    Omega2(m,n) = (At{m+1,2}-At{m+1,3}*At{end-1,5})/(At{end-1,5}+At{end,2});
                end
                Betas(1:length(Ap(:,n)),n) = Amdl.coeffs([2 4 6 8 10 12 16]);

                R2(n) = 1 - (At{end-1,2}/At{end,2});
            end
            ne = ne + 1 ;
            res_anova(ne).Ap = Ap;
            res_anova(ne).PEVs = PEVs;
            res_anova(ne).Fs= Fs;
            res_anova(ne).Omega2 = Omega2;
            res_anova(ne).Betas = Betas;
            res_anova(ne).R2 = R2;
            
            res_anova(ne).neurons_area = SPK.neurons_area ;
            res_anova(ne).neurons_info = SPK.neurons_info ;
            res_anova(ne).neurons_info_header = SPK.neurons_info_header ;
            res_anova(ne).session = list(s).name(1:8);   
            res_anova(ne).pref_ch = pref_ch(takebhv,:);
            
        end
    end
end
save([path2go 'res_ANOVAprefbin_' monkey '_v1.mat'],'res_anova','list','-v7.3')
%save([path2go 'res_ANOVAprefbin_REF_' monkey '_v1.mat'],'res_anova','list','-v7.3')


%%

clear
%path2go = '/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-BehavPrefChange/data/neurons/subset/'; %- path where SPKpool files are!
path2go = ('C:\Users\Fred\Dropbox\Rudebeck Lab\ANA-POTT-BehavPrefChange\data\neurons\subset\')

load([path2go 'res_ANOVAprefbin_X_v1.mat'])
%load([path2go 'res_ANOVAprefbin_M_v1.mat'])
%load([path2go 'res_ANOVAprefbin_REF_X_v1.mat'])
%load([path2go 'res_ANOVAprefbin_REF_M_v1.mat'])

area2test = {'vlPFC' 'OFC' 'IFG' 'LAI' 'Striatum' 'AMG'};

all_diff = [];
all_diff_omega = [];
all_diff_betas = [];
all_diff_PEVs = [];
all_sess = [];
all_units = [];
all_mk = [];
all_pref = [];
all_pref_center = [];
for s = 1 : length(res_anova)
    if ~isempty(res_anova(s).neurons_area)
        all_diff_omega = [all_diff_omega , res_anova(s).Omega2];
        all_diff_betas = [all_diff_betas , res_anova(s).Betas];
        all_diff_PEVs = [all_diff_PEVs , res_anova(s).PEVs];
        all_diff = [all_diff , res_anova(s).Ap];
        all_units = [all_units ; res_anova(s).neurons_area];
        all_mk = [all_mk ; repmat({res_anova(s).session(1)},size(res_anova(s).neurons_area,1),1)];
        all_sess = [all_sess ; repmat(s,size(res_anova(s).neurons_area,1),1)];
  %      all_pref = [all_pref ; repmat(res_anova(s).p_trend,size(res_anova(s).neurons_area,1),1)];
  %      all_pref_center = [all_pref_center ; repmat(res_anova(s).pref_center,size(res_anova(s).neurons_area,1),1)];
    end
end

keep = (sum(isnan(all_diff))==0)';


params = [3 1];
area_list = utils_POTT_areas;

for ar = 1 : length(area2test)
    eval(['takeit = ismember(all_units,area_list.' area2test{ar} ') & keep ;'])
    sigParam(:,ar) = sum(all_diff(:,takeit)<0.01,2);
    totalParam(1,ar) = length(all_diff(:,takeit));
end
figure;imagesc(sigParam./repmat(totalParam,size(sigParam,1),1))
          

% out of the proba neurons
clear sigParam
% params = {'Ju' 'Satiety' 'ChPb x Ju' 'ChPb x Sat' ' Ju x Sat' 'ChPb x Ju x Sat'}
params = {'ChPb' 'Satiety' 'ChPb x Ju' 'ChPb x Sat' ' Ju x Sat' 'ChPb x Ju x Sat'}
for ar = 1 : length(area2test)
    eval(['takeit = ismember(all_units,area_list.' area2test{ar} ') & keep ;'])
    probaUnits = sum(all_diff([1 4 5 7],:)<0.01)~=0 ;
    sigParam(:,ar) = sum(all_diff(2:end,takeit & probaUnits')<0.01,2);
    totalParam(1,ar) = length(all_diff(:,takeit & probaUnits'));
  %  juiceUnits = sum(all_diff([2 4 6 7],:)<0.01)~=0 ;
  %  sigParam(:,ar) = sum(all_diff([1 3:end],takeit & juiceUnits')<0.01,2);
  %  totalParam(1,ar) = length(all_diff(:,takeit & juiceUnits'));
    area2test_nb{ar} = [area2test{ar} '-' num2str(totalParam(1,ar))];
end
mat2plot = 100*(sigParam./repmat(totalParam,size(sigParam,1),1));
figure;imagesc(mat2plot)
for i = 1 : size(mat2plot,1)
    for j = 1 : size(mat2plot,2)
        text(j,i,num2str(mat2plot(i,j)))
    end
end
    set(gca,'Ytick', 1 : length(params),'YtickLabel',  params, ...
        'Xtick', 1 : length(area2test),'XtickLabel',  area2test_nb)
    





%% Pref ANOVA on SAME JUICE TRIALS - neurons


clear

monkey = 'X'
norm_me = @(data) -1+((data-min(data))*2)/(max(data)-min(data)) ;
norm_pb = @(data) -1+((data-0.1)*2)/(0.9-0.1) ;
predic_INS={'I_chosenproba'; 'I_unchosenproba' ;'I_chosenjuice'}; %- select any param you want to test.. order matters!

% path2go = '/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-BehavPrefChange/data/';
% %path2go = 'C:\\Users\fred\Dropbox\Rudebeck Lab\ANA-POTT-BehavPrefChange\data\';
% load([path2go 'Mimic_behav_norm_ALL_prevJuice_only.mat'])

path2go = '/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-BehavPrefChange/data/neurons/subset/'; %- path where SPKpool files are!
path2go2 = '/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-BehavPrefChange/data/';
load([path2go2 'Pref_bins.mat'])

if strcmp(monkey,'X')
    load([path2go2 'Mimic_behav_bins.mat'])
 %   load([path2go2 'Mimic_behav_bins_perms.mat'],'ALL','param')
    list = dir([path2go 'X*a_SPKpool.mat']);
    pref_ch = pref_ch_mk{2};
else
    load([path2go2 'Morbier_behav_bins.mat'])
 %   load([path2go2 'Morbier_behav_bins_perms.mat'],'ALL','param')
    list = dir([path2go 'M*a_SPKpool.mat']);
    pref_ch = pref_ch_mk{1};
end

clear days
for i = 1 : length(ALL)
    days{i,:} = [ALL(i).name(end-17:end-11) 'a'];
end

x = 0;

rmv = 100; % remove xxx ms on each side on every events (avoid overlaps between bins and smoothing problems)
% times_evts = {'FixFP_onset' 'Stim_onset' 'Resp_onset' 'FixResp_onset' 'FB_onset' 'Rew_onset' 'FB_offset'};
bins4decoding=[2]; %- perform decoding on subset on bins (stim and go period here)
binwindow = [200 800];
area_list = utils_POTT_areas;
area2test = {'a11ml' 'a13ml' 'a12mr' 'a12ol' 'IFG' 'LAI' 'AMG' 'Cd' 'PUT' }; % change order
area2test = {'vlPFC' 'OFC' 'IFG' 'LAI' 'Striatum' 'AMG'};

minNeurons = 5;

all_unit = cell(length(area2test),length(list),3);
all_unit_ft = cell(length(area2test),length(list),3);
all_unit_rt = cell(length(area2test),length(list),3);
pref_sig = NaN(length(list),4);
pca_explained = cell(length(area2test),length(list),3);
all_unit_data=cell(length(area2test),length(list),3);

meannorm  = @(data) (data-mean(mean(data)))/(max(max(data))-min(min(data)));
minnorm  = @(data) (data-min(min(data)))/(max(max(data))-min(min(data)));
normal = true;

final_data=table();

measures = {'all'} ;
clear res_anova 
ne = 0;
for s =  1 : length(list)
    disp(s)
    SPK = load([path2go list(s).name]);
    idxSess = ismember(days,list(s).name(1:8));
    takebhv = ismember(pref_ch(:,1),list(s).name(1:8));
    if ~isempty(SPK.neurons_area) && sum(takebhv)==1  %- skip sessions where no recorded neurons where in area2test!
        % spkfile = [ ALL(s).name(end-17:end-10) '_SPKpool.mat']
        % areafile = [ ALL(s).name(end-17:end-10) '_AreaInfo.mat']
        
        %load([path2go '/neurons/M072718a_SPKpool.mat'])
        
        %- find times
        bins2remove = (rmv/SPK.subsp)-1;
        
        n_evt = length(SPK.times_evts);
        lin=length(SPK.neurons_rank)*n_evt;
        nTimes = [0 ; (sum(abs(SPK.times(:,:)),2)/SPK.subsp)];
        col = sum(nTimes);
        bins = NaN(1,col);
        time = NaN(1,col);
        
        for b = 1 : length(SPK.times_evts)
            time(sum(nTimes(1:b))+1:sum(nTimes(1:b+1))) =  (SPK.times(b,1):SPK.subsp:SPK.times(b,2)-(SPK.subsp/2));
            bins(sum(nTimes(1:b))+1:sum(nTimes(1:b+1))) = b;
            bins([sum(nTimes(1:b))+1:sum(nTimes(1:b))+1+bins2remove   sum(nTimes(1:b+1))-bins2remove:sum(nTimes(1:b+1))    ]) = NaN; %- remove 100 ms each side (avoid overlaps between bins and smoothing problems)
        end
        
        %- keep only considered bins
        % SPK.SPK_INS(:,~ismember(bins,bins4decoding))=[];
        time2cons = ismember(bins,bins4decoding) & time>=binwindow(1) & time<=binwindow(2);
        SPK.SPK_INS = SPK.SPK_INS(:,time2cons);
        time(~time2cons) = [];
        
        
        %  SPK.SPK_INS = SPK.SPK_INS(:,296:316);
        
        
        %- find trial ID included for each 30-trial behavioral model
%         clear trID
%         dumm = find(ALL(idxSess).trial2take);
%         for b = 1 : length(dumm)-param.binSize
%             trID(b,:) = dumm(b:b+param.binSize);
%         end
        
        % find trial ID for the neurons
        trID_neurons = find(ALL(idxSess).TrialType(1,:)==2 & ALL(idxSess).TrialType(2,:)==0);
        
        %     figure;plot(ALL(s).TrialType(11,trID_neurons));hold on
        %     plot(SPK.TrialType_INS(12,:))
        [ALL(idxSess).Z_trend_ft,~]=Mann_Kendall(ALL(idxSess).ft_bins,0.01);
        [ALL(idxSess).Z_trend_rt,~]=Mann_Kendall(ALL(idxSess).rt_bins,0.01);
        [ALL(idxSess).Z_trend_pref,~]=Mann_Kendall(ALL(idxSess).pref_bins,0.01);
        
        if length(trID_neurons)~=size(SPK.TrialType_INS,2)
            warning('Problem in # of trials......')
            x = x + 1;
            pb_files(x)=s;

        else

            pref_ch(takebhv,:)

            dumm = SPK.TrialType_INS(ismember(SPK.TrialType_header,{'nTr' 'I_chosenproba' 'I_chosenjuice' }),:)';
            dumm = [dumm(:,1) dumm(:,2)+(100*dumm(:,3))];

            when = @(arg,cd) SPK.TrialType_INS(strcmp(SPK.TrialType_header,arg),:)==cd ;
            diff_juice = (when('I_juiceL',1) & when('I_juiceR',2)) | (when('I_juiceL',2) & when('I_juiceR',1)); %- take only the different juice trials
            same_juice = (when('I_juiceL',1) & when('I_juiceR',1)) | (when('I_juiceL',2) & when('I_juiceR',2)) ; %- take same juice trials

            %- extract only the diff juice trials
            tt = find(diff_juice);

            keep = [130 150 170 190 230 250 270 290];
            %  remove = ~ismember(dumm(:,2),keep) | ~diff_juice';

            %- pre period
            preTr = pref_ch{takebhv,2}; %- out of the diff juice trials
            preTr4spk = tt(preTr); %- trial to take for neurons

            % sum(diff( pref_ch{takebhv,3})~=1) %- used to check if sometimes gaps between multiple post periods.. Not the case!
            %- because no gaps, can do simple boundaries around min max tr nb
            pre_boundaries = [min(preTr4spk) max(preTr4spk)];
            pre_period = trID_neurons>=trID_neurons(pre_boundaries(1)) & trID_neurons<=trID_neurons(pre_boundaries(2));
            preTr_same = find(same_juice & pre_period);

            preTr_fact = dumm(preTr_same,:);
            preRemove = ~ismember(preTr_fact(:,2),keep) ;
            preTr_fact(preRemove,:)=[];

            %- post period
            postTr = pref_ch{takebhv,3}; %- out of the diff juice trials
            postTr4spk = tt(postTr); %- trial to take for neurons

            % sum(diff( pref_ch{takebhv,3})~=1) %- used to check if sometimes gaps between multiple post periods.. Not the case!
            %- because no gaps, can do simple boundaries around min max tr nb
            post_boundaries = [min(postTr4spk) max(postTr4spk)];
            post_period = trID_neurons>=trID_neurons(post_boundaries(1)) & trID_neurons<=trID_neurons(post_boundaries(2));
            postTr_same = find(same_juice & post_period);

            postTr_fact = dumm(postTr_same,:);
            postRemove = ~ismember(postTr_fact(:,2),keep) ;
            postTr_fact(postRemove,:)=[];

            if pref_ch{takebhv,4}(2)-pref_ch{takebhv,4}(1)>0
                allfactors = [preTr_fact(:,1) norm_pb(mod(preTr_fact(:,2),100)/100) abs(floor(preTr_fact(:,2)/100)-3) ones(size(preTr_fact(:,1))) ; ...
                                postTr_fact(:,1) norm_pb(mod(postTr_fact(:,2),100)/100) abs(floor(postTr_fact(:,2)/100)-3) 2*ones(size(postTr_fact(:,1)))]; %- inverse juice (1 = valued / 2 = devalued)
            else
                allfactors = [preTr_fact(:,1) norm_pb(mod(preTr_fact(:,2),100)/100) floor(preTr_fact(:,2)/100) ones(size(preTr_fact(:,1))) ; ...
                              postTr_fact(:,1) norm_pb(mod(postTr_fact(:,2),100)/100) floor(postTr_fact(:,2)/100) 2*ones(size(postTr_fact(:,1)))]; %- no inverse juice (1 = valued / 2 = devalued)
            end
            

            clear Ap PEVs Fs Omega2 R2 Betas
            for n = 1 : length(SPK.neurons_area)
                
                data = mean(SPK.SPK_INS(SPK.Tr_Clust_INS(:,2)==n,:),2);
                
                [Ap(:,n),At,Amdl,~] = anovan(data(allfactors(:,1)),allfactors(:,2:4),'continuous',[1],...
                    'varname',{'ChosenProba' 'Val_Deval' 'Pre_Post'},...
                    'model',[1 0 0;0 1 0;0 0 1; 1 1 0 ; 1 0 1 ; 0 1 1 ; 1 1 1],'display','off');

                for m = 1 : length(Ap(:,n))
                    PEVs(m,n) = ((At{m+1,2}/At{end,2}))*100;
                    Fs(m,n) = At{m+1,6};
                    Omega2(m,n) = (At{m+1,2}-At{m+1,3}*At{end-1,5})/(At{end-1,5}+At{end,2});
                end
                Betas(1:length(Ap(:,n)),n) = Amdl.coeffs([2 4 6 8 10 12 16]);

                R2(n) = 1 - (At{end-1,2}/At{end,2});
            end
            ne = ne + 1 ;
            res_anova(ne).Ap = Ap;
            res_anova(ne).PEVs = PEVs;
            res_anova(ne).Fs= Fs;
            res_anova(ne).Omega2 = Omega2;
            res_anova(ne).Betas = Betas;
            res_anova(ne).R2 = R2;
            
            res_anova(ne).neurons_area = SPK.neurons_area ;
            res_anova(ne).neurons_info = SPK.neurons_info ;
            res_anova(ne).neurons_info_header = SPK.neurons_info_header ;
            res_anova(ne).session = list(s).name(1:8);   
            res_anova(ne).pref_ch = pref_ch(takebhv,:);
            
        end
    end
end
save([path2go 'res_ANOVAprefbin_same_' monkey '_v1.mat'],'res_anova','list','-v7.3')
%save([path2go 'res_ANOVAprefbin_REF_' monkey '_v1.mat'],'res_anova','list','-v7.3')

   
          

%%



clear
%path2go = '/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-BehavPrefChange/data/neurons/subset/'; %- path where SPKpool files are!
path2go = ('C:\Users\Fred\Dropbox\Rudebeck Lab\ANA-POTT-BehavPrefChange\data\neurons\subset\')

X = load([path2go 'res_ANOVAprefbin_X_v1.mat'])
M = load([path2go 'res_ANOVAprefbin_M_v1.mat'])
%load([path2go 'res_ANOVAprefbin_REF_X_v1.mat'])
%load([path2go 'res_ANOVAprefbin_REF_M_v1.mat'])
res_anova = [X.res_anova M.res_anova]
area2test = {'vlPFC' 'OFC' 'IFG' 'LAI' 'Striatum' 'AMG'};

all_diff = [];
all_diff_omega = [];
all_diff_betas = [];
all_diff_PEVs = [];
all_sess = [];
all_units = [];
all_mk = [];
all_pref = [];
all_pref_center = [];
for s = 1 : length(res_anova)
    if ~isempty(res_anova(s).neurons_area)
        all_diff_omega = [all_diff_omega , res_anova(s).Omega2];
        all_diff_betas = [all_diff_betas , res_anova(s).Betas];
        all_diff_PEVs = [all_diff_PEVs , res_anova(s).PEVs];
        all_diff = [all_diff , res_anova(s).Ap];
        all_units = [all_units ; res_anova(s).neurons_area];
        all_mk = [all_mk ; repmat({res_anova(s).session(1)},size(res_anova(s).neurons_area,1),1)];
        all_sess = [all_sess ; repmat(s,size(res_anova(s).neurons_area,1),1)];
  %      all_pref = [all_pref ; repmat(res_anova(s).p_trend,size(res_anova(s).neurons_area,1),1)];
  %      all_pref_center = [all_pref_center ; repmat(res_anova(s).pref_center,size(res_anova(s).neurons_area,1),1)];
    end
end

keep = (sum(isnan(all_diff))==0)';


params = [3 1];
area_list = utils_POTT_areas;

for ar = 1 : length(area2test)
    eval(['takeit = ismember(all_units,area_list.' area2test{ar} ') & keep ;'])
    sigParam(:,ar) = sum(all_diff(:,takeit)<0.01,2);
    totalParam(1,ar) = length(all_diff(:,takeit));
end
figure;imagesc(sigParam./repmat(totalParam,size(sigParam,1),1))
          

% out of the proba neurons
clear sigParam
% params = {'Ju' 'Satiety' 'ChPb x Ju' 'ChPb x Sat' ' Ju x Sat' 'ChPb x Ju x Sat'}
 params = {'ChPb' 'Satiety' 'ChPb x Ju' 'ChPb x Sat' ' Ju x Sat' 'ChPb x Ju x Sat'}
for ar = 1 : length(area2test)
    eval(['takeit = ismember(all_units,area_list.' area2test{ar} ') & keep ;'])
  %  probaUnits = sum(all_diff([1 4 5 7],:)<0.01)~=0 ;
  %  sigParam(:,ar) = sum(all_diff(2:end,takeit & probaUnits')<0.01,2);
  %  totalParam(1,ar) = length(all_diff(:,takeit & probaUnits'));
    juiceUnits = sum(all_diff([2 4 6 7],:)<0.01)~=0 ;
    sigParam(:,ar) = sum(all_diff([1 3:end],takeit & juiceUnits')<0.01,2);
    totalParam(1,ar) = length(all_diff(:,takeit & juiceUnits'));
    area2test_nb{ar} = [area2test{ar} '-' num2str(totalParam(1,ar))];
end
mat2plot = 100*(sigParam./repmat(totalParam,size(sigParam,1),1));
figure;imagesc(mat2plot)
for i = 1 : size(mat2plot,1)
    for j = 1 : size(mat2plot,2)
        text(j,i,num2str(mat2plot(i,j)))
    end
end
    set(gca,'Ytick', 1 : length(params),'YtickLabel',  params, ...
        'Xtick', 1 : length(area2test),'XtickLabel',  area2test_nb)
    

% out of the proba neurons
clear sigSat_proba sigSat_juice sigSat_all
% params = {'Ju' 'Satiety' 'ChPb x Ju' 'ChPb x Sat' ' Ju x Sat' 'ChPb x Ju x Sat'}
 params = {'ChPb' 'Satiety' 'ChPb x Ju' 'ChPb x Sat' ' Ju x Sat' 'ChPb x Ju x Sat'}
for ar = 1 : length(area2test)
    eval(['takeit = ismember(all_units,area_list.' area2test{ar} ') & keep ;'])

    probaUnits = sum(all_diff([1 4 5 7],:)<0.01)~=0 ;
   % probaUnits = (all_diff([1 ],:)<0.01) ;
    sigSat_proba(ar,:) = [sum(sum(all_diff([3 5 6 7],takeit & probaUnits')<0.01)~=0)  length(all_diff(:,takeit & probaUnits'))];
    sigSat_proba_M(ar,:) = [sum(sum(all_diff([3 5 6 7],takeit & probaUnits' & ismember(all_mk,'M'))<0.01)~=0)  length(all_diff(:,takeit & probaUnits' & ismember(all_mk,'M')))];
    sigSat_proba_X(ar,:) = [sum(sum(all_diff([3 5 6 7],takeit & probaUnits' & ismember(all_mk,'X'))<0.01)~=0)  length(all_diff(:,takeit & probaUnits' & ismember(all_mk,'X')))];

    juiceUnits = sum(all_diff([2 4 6 7],:)<0.01)~=0 ;
   % juiceUnits = (all_diff([2],:)<0.01) ;
    sigSat_juice(ar,:) = [sum(sum(all_diff([3 5 6 7],takeit & juiceUnits')<0.01)~=0)  length(all_diff(:,takeit & juiceUnits'))];
    sigSat_juice_M(ar,:) = [sum(sum(all_diff([3 5 6 7],takeit & juiceUnits' & ismember(all_mk,'M'))<0.01)~=0)  length(all_diff(:,takeit & juiceUnits' & ismember(all_mk,'M')))];
    sigSat_juice_X(ar,:) = [sum(sum(all_diff([3 5 6 7],takeit & juiceUnits' & ismember(all_mk,'X'))<0.01)~=0)  length(all_diff(:,takeit & juiceUnits' & ismember(all_mk,'X')))];

    smthgUnits = sum(all_diff([1 2 3 4 5 6 7],:)<0.01)~=0 ;
    sigSat_smthg(ar,:) = [sum(sum(all_diff([3 5 6 7],takeit & smthgUnits')<0.01)~=0)  length(all_diff(:,takeit & smthgUnits'))];
    sigSat_smthg_M(ar,:) = [sum(sum(all_diff([3 5 6 7],takeit & smthgUnits' & ismember(all_mk,'M'))<0.01)~=0)  length(all_diff(:,takeit & smthgUnits' & ismember(all_mk,'M')))];
    sigSat_smthg_X(ar,:) = [sum(sum(all_diff([3 5 6 7],takeit & smthgUnits' & ismember(all_mk,'X'))<0.01)~=0)  length(all_diff(:,takeit & smthgUnits' & ismember(all_mk,'X')))];


    sigSat_all(ar,:) = [sum(sum(all_diff([3 5 6 7],takeit)<0.01)~=0)  length(all_diff(:,takeit ))];
end
figure;bar([sigSat_proba(:,1)./sigSat_proba(:,2) ...
            sigSat_juice(:,1)./sigSat_juice(:,2) ...
            sigSat_smthg(:,1)./sigSat_smthg(:,2)  ]); hold on
for i = 1 : length(area2test)
    text(i-.225,sigSat_proba(i,1)./sigSat_proba(i,2)+.01,[num2str(sigSat_proba(i,1)) '/' num2str(sigSat_proba(i,2))],'Rotation',90)
    text(i,sigSat_juice(i,1)./sigSat_juice(i,2)+.01,[num2str(sigSat_juice(i,1)) '/' num2str(sigSat_juice(i,2))],'Rotation',90)
    text(i+.225,sigSat_smthg(i,1)./sigSat_smthg(i,2)+.01,[num2str(sigSat_smthg(i,1)) '/' num2str(sigSat_smthg(i,2))],'Rotation',90)
end
    set(gca,'Xtick', 1 : length(area2test),'XtickLabel',  area2test)
    

sigSat_proba(:,1)./sigSat_proba(:,2) 
sigSat_juice(:,1)./sigSat_juice(:,2) 
sigSat_all(:,1)./sigSat_all(:,2) 




sigSat_juice_M


        %- binomial model for perc Sig Neurons 
        x1= [];x2=[];x3 = [];
        for mm = 1 : 2
            if mm ==1
                nN = sigSat_proba_M;
             %   nN = sigSat_juice_M;
                nN = sigSat_smthg_M;
            else
                nN = sigSat_proba_X;
             %   nN = sigSat_juice_X;
                nN = sigSat_smthg_X;
            end
            for i  = 1 : length(nN(:,1))
                x1 = [x1 ; repmat(area2test(i),nN(i,2),1)];
            end
            for i  = 1 : length(nN(:,1))
                x2 = [x2 ; repmat(true,nN(i,1),1) ; repmat(false,nN(i,2)-nN(i,1),1)];
            end
            x3 = [x3 ; repmat(['mk' num2str(mm)],sum(nN(:,2)),1)];
        end
        
        modeldata = table(logical(x2),x1,x3,'VariableNames',{'sig' 'area' 'mk'})
        
        lme_percSig = fitglme(modeldata,'sig ~ 1 + area + (1|mk) ','Distribution','binomial','Link','logit')
        
        anova(lme_percSig)
        [~,wald_percNeurons,~,pval_adj_percNeurons] = area_posthoc(lme_percSig,area2test,'y');












