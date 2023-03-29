%% SPK ANOVA for preference modulation
%- Require the time windows for PRE/POST pref change extracted in POTT_BHV_001_GLM

%- Author: Fred M. Stoll, Icahn School of Medicine at Mount Sinai, NY
%- Date: 2022.12

%% Pref ANOVA - neurons

clear

monkey = 'X';
norm_me = @(data) -1+((data-min(data))*2)/(max(data)-min(data)) ;
norm_pb = @(data) -1+((data-0.1)*2)/(0.9-0.1) ;
predic_INS={'I_chosenproba'; 'I_unchosenproba' ;'I_chosenjuice'}; %- select any param you want to test.. order matters!

% path2go = '/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-BehavPrefChange/data/';
path2go = 'C:\\Users\fred\Dropbox\Rudebeck Lab\ANA-POTT-BehavPrefChange\data\';
% load([path2go 'Mimic_behav_norm_ALL_prevJuice_only.mat'])

%path2go = '/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-BehavPrefChange/data/neurons/subset/'; %- path where SPKpool files are!
%path2go2 = '/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-BehavPrefChange/data/';
path2go = 'C:\\Users\fred\Dropbox\Rudebeck Lab\ANA-POTT-BehavPrefChange\data\neurons\subset\';
path2go2 = 'C:\\Users\fred\Dropbox\Rudebeck Lab\ANA-POTT-BehavPrefChange\data\';
%load([path2go2 'Pref_bins.mat'])

load([path2go2 'FT_bins.mat'])
%- remove sessions where PREF also change..
PREF = load([path2go2 'Pref_bins.mat']);
pref_ch_mk{1}(ismember(pref_ch_mk{1}(:,1),PREF.pref_ch_mk{1}(:,1)),:)=[];
pref_ch_mk{2}(ismember(pref_ch_mk{2}(:,1),PREF.pref_ch_mk{2}(:,1)),:)=[];

%sum(ismember(PREF.pref_ch_mk{1}(:,1),pref_ch_mk{1}(:,1))) 
%length(ismember(PREF.pref_ch_mk{1}(:,1),pref_ch_mk{1}(:,1))) 

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
nbPerm = 100;

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
            

            clear Ap PEVs Fs Omega2 R2 Betas FR factors PEVs_perm Omega2_perm Ap_perm R2_perm

            for n = 1 : length(SPK.neurons_area)
                
                data = mean(SPK.SPK_INS(SPK.Tr_Clust_INS(:,2)==n,:),2);
                
                [Ap(:,n),At,Amdl,~] = anovan(data(allfactors(:,1)),allfactors(:,2:4),'continuous',[1],'varname',{'ChosenProba' 'Val_Deval' 'Pre_Post'},'model',[1 0 0;0 1 0;0 0 1; 1 1 0 ; 1 0 1 ; 0 1 1 ; 1 1 1],'display','off');
                par = 1 : length(Ap(:,n));
                PEVs(:,n) = (([At{par+1,2}]/At{end,2}))*100;
                Fs(:,n) = [At{par+1,6}];
                Omega2(:,n) = ([At{par+1,2}]-[At{par+1,3}]*At{end-1,5})/(At{end-1,5}+At{end,2});
                Betas(1:length(Ap(:,n)),n) = Amdl.coeffs([2 4 6 8 10 12 16]);

                R2(n) = 1 - (At{end-1,2}/At{end,2});

                FR(n,:) = data(allfactors(:,1));
                factors = allfactors(:,2:4);

                %- permutations
                parfor pp = 1 : nbPerm
                    ordr = randperm(length(allfactors(:,1)));
                    [Ap_perm(:,n,pp),At_perm,~,~] = anovan(FR(n,ordr),factors,'continuous',[1],'varname',{'ChosenProba' 'Val_Deval' 'Pre_Post'},'model',[1 0 0;0 1 0;0 0 1; 1 1 0 ; 1 0 1 ; 0 1 1 ; 1 1 1],'display','off');
                    par = 1 : length(Ap(:,n));
                    PEVs_perm(:,n,pp) = (([At_perm{par+1,2}]/At_perm{end,2}))*100;
                    Omega2_perm(:,n,pp) = ([At_perm{par+1,2}]-[At_perm{par+1,3}]*At_perm{end-1,5})/(At_perm{end-1,5}+At_perm{end,2});
                    R2_perm(n,pp) = 1 - (At_perm{end-1,2}/At_perm{end,2});
                end
                
            end
            
            ne = ne + 1 ;
            res_anova(ne).Ap = Ap;
            res_anova(ne).Ap_perm = Ap_perm;
            res_anova(ne).PEVs = PEVs;
            res_anova(ne).PEVs_perm = PEVs_perm;
            res_anova(ne).Fs= Fs;
            res_anova(ne).Omega2 = Omega2;
            res_anova(ne).Omega2_perm = Omega2_perm;
            res_anova(ne).Betas = Betas;
            res_anova(ne).R2 = R2;
            res_anova(ne).R2_perm = R2_perm;
            res_anova(ne).FR = FR;
            res_anova(ne).factors = factors;
            
            res_anova(ne).neurons_area = SPK.neurons_area ;
            res_anova(ne).neurons_info = SPK.neurons_info ;
            res_anova(ne).neurons_info_header = SPK.neurons_info_header ;
            res_anova(ne).session = list(s).name(1:8);   
            res_anova(ne).pref_ch = pref_ch(takebhv,:);
            
        end
    end
end
%save([path2go 'res_ANOVAprefbin_' monkey '_v1.mat'],'res_anova','list','-v7.3')
save([path2go 'res_ANOVAftbin_' monkey '_v1.mat'],'res_anova','list','-v7.3')
%save([path2go 'res_ANOVAprefbin_REF_' monkey '_v1.mat'],'res_anova','list','-v7.3')


%% Post hoc

clear
path2go = '/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-BehavPrefChange/data/neurons/subset/'; %- path where SPKpool files are!
%path2go = ('C:\Users\Fred\Dropbox\Rudebeck Lab\ANA-POTT-BehavPrefChange\data\neurons\subset\')

sigAnova = load('C:\Users\Fred\Dropbox\Rudebeck Lab\ANA-POTT-BehavPrefChange\POTT_sigUnits_name.mat')

X = load([path2go 'res_ANOVAprefbin_X_v1.mat'])
M = load([path2go 'res_ANOVAprefbin_M_v1.mat'])
%X = load([path2go 'res_ANOVAftbin_X_v1.mat'])
%M = load([path2go 'res_ANOVAftbin_M_v1.mat'])
%load([path2go 'res_ANOVAprefbin_REF_X_v1.mat'])
%load([path2go 'res_ANOVAprefbin_REF_M_v1.mat'])
res_anova = [X.res_anova M.res_anova];
area2test = {'vlPFC' 'OFC' 'IFG' 'LAI' 'Striatum' 'AMG'};
area_list = utils_POTT_areas;

thr_sig = 0.01;

%- combine all sessions
all_diff = [res_anova(:).Ap];
all_diff_perm = [res_anova(:).Ap_perm];
all_units = cat(1,res_anova(:).neurons_area);
all_beta = [res_anova(:).Betas];

all_mk = [];
for s = 1 : length(res_anova)
    all_mk = [all_mk ; repmat({res_anova(s).session(1)},size(res_anova(s).neurons_area,1),1)];
end

%- remove neurons with ANOVAs not converging
keep = (sum(isnan(all_diff))==0)';

all_names = [];
for s = 1 : length(res_anova)
    for n = 1 : length(res_anova(s).neurons_area)
    if res_anova(s).neurons_info(n,1)<10; addzeros = '00';
    elseif res_anova(s).neurons_info(n,1)<100; addzeros = '0';
    else addzeros = '';
    end
    if res_anova(s).neurons_info(n,2)<10; addzeros2 = '0';
    else addzeros2 = '';
    end

    all_names = [all_names ; [res_anova(s).session '_Ch' addzeros num2str(res_anova(s).neurons_info(n,1)) ...
        '_Clus' addzeros2 num2str(res_anova(s).neurons_info(n,2)) ]];
    end
end

%- find which neurons were sig for proba or juice in main ana...
probaUnits = ismember(cellstr(all_names), cellstr(sigAnova.name_sig{2}))';
juiceUnits = ismember(cellstr(all_names), cellstr(sigAnova.name_sig{1}))';

% for ar = 1 : length(area2test)
%     eval(['takeit = ismember(all_units,area_list.' area2test{ar} ') & keep ;'])
%     sigParam(:,ar) = sum(all_diff(:,takeit)<0.01,2);
%     totalParam(1,ar) = length(all_diff(:,takeit));
% end
% figure;imagesc(sigParam./repmat(totalParam,size(sigParam,1),1))
          

% % out of the proba neurons
% clear sigParam
% % params = {'Ju' 'Satiety' 'ChPb x Ju' 'ChPb x Sat' ' Ju x Sat' 'ChPb x Ju x Sat'}
%  params = {'ChPb' 'Satiety' 'ChPb x Ju' 'ChPb x Sat' ' Ju x Sat' 'ChPb x Ju x Sat'}
% for ar = 1 : length(area2test)
%     eval(['takeit = ismember(all_units,area_list.' area2test{ar} ') & keep ;'])
%   %  probaUnits = sum(all_diff([1 4 5 7],:)<0.01)~=0 ;
%   %  sigParam(:,ar) = sum(all_diff(2:end,takeit & probaUnits')<0.01,2);
%   %  totalParam(1,ar) = length(all_diff(:,takeit & probaUnits'));
%     juiceUnits = sum(all_diff([2 4 6 7],:)<0.01)~=0 ;
%     sigParam(:,ar) = sum(all_diff([1 3:end],takeit & juiceUnits')<0.01,2);
%     totalParam(1,ar) = length(all_diff(:,takeit & juiceUnits'));
%     area2test_nb{ar} = [area2test{ar} '-' num2str(totalParam(1,ar))];
% end
% mat2plot = 100*(sigParam./repmat(totalParam,size(sigParam,1),1));
% figure;imagesc(mat2plot)
% for i = 1 : size(mat2plot,1)
%     for j = 1 : size(mat2plot,2)
%         text(j,i,num2str(mat2plot(i,j)))
%     end
% end
%     set(gca,'Ytick', 1 : length(params),'YtickLabel',  params, ...
%         'Xtick', 1 : length(area2test),'XtickLabel',  area2test_nb)
    

% out of the proba/juice/smthg/all neurons
clear sigSat_proba sigSat_juice sigSat_all satUnits satUnits_M satUnits_X prop_SatIncr prop_SatIncr_X prop_SatIncr_M
mat_ref = all_diff; % use the true 'proba', 'juice' etc neurons to see proportion sig in permutations
for ar = 1 : length(area2test)
    eval(['takeit = ismember(all_units,area_list.' area2test{ar} ') & keep ;'])
    
    for pp = 1 : size(all_diff_perm,3)+1
        if pp == 1
            mat = all_diff;
        else
            mat = squeeze(all_diff_perm(:,:,pp-1));
        end
        
     %   probaUnits = sum(mat_ref([1 4 5 7],:)<thr_sig)~=0 ;
        % probaUnits = (mat([1 ],:)<thr_sig) ;
        sigSat_proba(ar,:,pp) = [sum(sum(mat([3 5 6 7],takeit & probaUnits')<thr_sig)~=0)  length(mat(1,takeit & probaUnits'))];
        sigSat_proba_M(ar,:,pp) = [sum(sum(mat([3 5 6 7],takeit & probaUnits' & ismember(all_mk,'M'))<thr_sig)~=0)  length(mat(1,takeit & probaUnits' & ismember(all_mk,'M')))];
        sigSat_proba_X(ar,:,pp) = [sum(sum(mat([3 5 6 7],takeit & probaUnits' & ismember(all_mk,'X'))<thr_sig)~=0)  length(mat(1,takeit & probaUnits' & ismember(all_mk,'X')))];
        
    %    juiceUnits = sum(mat_ref([2 4 6 7],:)<thr_sig)~=0 ;
        % juiceUnits = (mat([2],:)<thr_sig) ;
        sigSat_juice(ar,:,pp) = [sum(sum(mat([3 5 6 7],takeit & juiceUnits')<thr_sig)~=0)  length(mat(1,takeit & juiceUnits'))];
        sigSat_juice_M(ar,:,pp) = [sum(sum(mat([3 5 6 7],takeit & juiceUnits' & ismember(all_mk,'M'))<thr_sig)~=0)  length(mat(1,takeit & juiceUnits' & ismember(all_mk,'M')))];
        sigSat_juice_X(ar,:,pp) = [sum(sum(mat([3 5 6 7],takeit & juiceUnits' & ismember(all_mk,'X'))<thr_sig)~=0)  length(mat(1,takeit & juiceUnits' & ismember(all_mk,'X')))];
        
      %  smthgUnits = sum(mat_ref([1 2 3 4 5 6 7],:)<thr_sig)~=0 ;
        smthgUnits = juiceUnits | probaUnits ;
        sigSat_smthg(ar,:,pp) = [sum(sum(mat([3 5 6 7],takeit & smthgUnits')<thr_sig)~=0)  length(mat(1,takeit & smthgUnits'))];
        sigSat_smthg_M(ar,:,pp) = [sum(sum(mat([3 5 6 7],takeit & smthgUnits' & ismember(all_mk,'M'))<thr_sig)~=0)  length(mat(1,takeit & smthgUnits' & ismember(all_mk,'M')))];
        sigSat_smthg_X(ar,:,pp) = [sum(sum(mat([3 5 6 7],takeit & smthgUnits' & ismember(all_mk,'X'))<thr_sig)~=0)  length(mat(1,takeit & smthgUnits' & ismember(all_mk,'X')))];
        
    if pp == 1
        sat_effect = sum(mat([3 5 6 7],:)<thr_sig)~=0;
        all_beta_sig = all_beta;
        all_beta_sig(mat>thr_sig)=NaN;
        dumm = sign(nanmean(all_beta_sig([3 5 6 7],takeit & smthgUnits' & sat_effect')));
        prop_SatIncr(ar,:) = [sum(dumm==1) length(dumm)];
        dumm = sign(nanmean(all_beta_sig([3 5 6 7],takeit & smthgUnits' & sat_effect' & ismember(all_mk,'M'))));
        prop_SatIncr_M(ar,:) = [sum(dumm==1) length(dumm)];
        dumm = sign(nanmean(all_beta_sig([3 5 6 7],takeit & smthgUnits' & sat_effect' & ismember(all_mk,'X'))));
        prop_SatIncr_X(ar,:) = [sum(dumm==1) length(dumm)];
    end
        nthgUnits = ~juiceUnits | ~probaUnits ;
        sigSat_nthg(ar,:,pp) = [sum(sum(mat([3 5 6 7],takeit & nthgUnits')<thr_sig)~=0)  length(mat(1,takeit & nthgUnits'))];
        sigSat_nthg_M(ar,:,pp) = [sum(sum(mat([3 5 6 7],takeit & nthgUnits' & ismember(all_mk,'M'))<thr_sig)~=0)  length(mat(1,takeit & nthgUnits' & ismember(all_mk,'M')))];
        sigSat_nthg_X(ar,:,pp) = [sum(sum(mat([3 5 6 7],takeit & nthgUnits' & ismember(all_mk,'X'))<thr_sig)~=0)  length(mat(1,takeit & nthgUnits' & ismember(all_mk,'X')))];

        sigSat_all(ar,:,pp) = [sum(sum(mat([3 5 6 7],takeit)<thr_sig)~=0)  length(mat(1,takeit ))];
        sigSat_all_M(ar,:,pp) = [sum(sum(mat([3 5 6 7],takeit & ismember(all_mk,'M'))<thr_sig)~=0)  length(mat(1,takeit & ismember(all_mk,'M')))];
        sigSat_all_X(ar,:,pp) = [sum(sum(mat([3 5 6 7],takeit & ismember(all_mk,'X'))<thr_sig)~=0)  length(mat(1,takeit & ismember(all_mk,'X')))];
        
        %- out of all the satiety neurons, what is the most prevalent pattern??
      %  satUnits(ar,:,pp) = [sum(mat_ref([3],takeit)<thr_sig & sum(mat_ref([5 6 7],takeit)<thr_sig)==0) ; ...
      %      sum(sum(mat_ref([5 7],takeit)<thr_sig)~=0) ; ...
      %      sum(sum(mat_ref([6 7],takeit)<thr_sig)~=0) ; ...
      %      sum(sum(mat_ref([3 5 6 7],takeit)<thr_sig)~=0)]';

        satUnits(ar,:,pp) = [sum(mat([3],takeit)<thr_sig & sum(mat([5 6 7],takeit)<thr_sig)==0) ; ...
            sum(mat([5],takeit)<thr_sig & sum(mat([6 7],takeit)<thr_sig)==0) ; ...
            sum(mat([6],takeit)<thr_sig & sum(mat([5 7],takeit)<thr_sig)==0) ; ...
            sum(mat([7],takeit)<thr_sig) ; ...
            sum(sum(mat([3 5 6 7],takeit)<thr_sig)~=0)]';
        
        satUnits_M(ar,:,pp) = [sum(mat([3],takeit & ismember(all_mk,'M'))<thr_sig & sum(mat([5 6 7],takeit & ismember(all_mk,'M'))<thr_sig)==0) ; ...
            sum(mat([5],takeit & ismember(all_mk,'M'))<thr_sig & sum(mat([6 7],takeit & ismember(all_mk,'M'))<thr_sig)==0) ; ...
            sum(mat([6],takeit & ismember(all_mk,'M'))<thr_sig & sum(mat([5 7],takeit & ismember(all_mk,'M'))<thr_sig)==0) ; ...
            sum(mat([7],takeit & ismember(all_mk,'M'))<thr_sig) ; ...
            sum(sum(mat([3 5 6 7],takeit & ismember(all_mk,'M'))<thr_sig)~=0)]';
        
        satUnits_X(ar,:,pp) = [sum(mat([3],takeit & ismember(all_mk,'X'))<thr_sig & sum(mat([5 6 7],takeit & ismember(all_mk,'X'))<thr_sig)==0) ; ...
            sum(mat([5],takeit & ismember(all_mk,'X'))<thr_sig & sum(mat([6 7],takeit & ismember(all_mk,'X'))<thr_sig)==0) ; ...
            sum(mat([6],takeit & ismember(all_mk,'X'))<thr_sig & sum(mat([5 7],takeit & ismember(all_mk,'X'))<thr_sig)==0) ; ...
            sum(mat([7],takeit & ismember(all_mk,'X'))<thr_sig) ; ...
            sum(sum(mat([3 5 6 7],takeit & ismember(all_mk,'X'))<thr_sig)~=0)]';
        
    end
end

%- check against perm.
test = {'proba' 'juice' 'smthg' 'nthg' 'all'};
for t = 1 : length(test)
    eval(['sigSat_true = sigSat_' test{t} '(:,1,1)./sigSat_' test{t} '(:,2,1);'])
    eval(['sigSat_perm = squeeze(sigSat_' test{t} '(:,1,2:end))./squeeze(sigSat_' test{t} '(:,2,2:end));'])
    pval(:,t) = sum(sigSat_true<=sigSat_perm,2)/size(sigSat_perm,2); %- one-sided test
end
% sigSat_perm = mean(squeeze(sigSat_proba(:,1,2:end))./squeeze(sigSat_proba(:,2,2:end)),2)

%- quick plot
% figure;bar([sigSat_proba(:,1,1)./sigSat_proba(:,2,1) ...
%             sigSat_juice(:,1,1)./sigSat_juice(:,2,1) ...
%             sigSat_smthg(:,1,1)./sigSat_smthg(:,2,1)  ]); hold on
% for i = 1 : length(area2test)
%     text(i-.225,sigSat_proba(i,1,1)./sigSat_proba(i,2,1)+.01,[num2str(sigSat_proba(i,1,1)) '/' num2str(sigSat_proba(i,2,1))],'Rotation',90)
%     text(i,sigSat_juice(i,1,1)./sigSat_juice(i,2,1)+.01,[num2str(sigSat_juice(i,1,1)) '/' num2str(sigSat_juice(i,2,1))],'Rotation',90)
%     text(i+.225,sigSat_smthg(i,1,1)./sigSat_smthg(i,2,1)+.01,[num2str(sigSat_smthg(i,1,1)) '/' num2str(sigSat_smthg(i,2,1))],'Rotation',90)
% end
% set(gca,'Xtick', 1 : length(area2test),'XtickLabel',  area2test)

%- COLOR ASSIGNMENT
order = [3 1 2 5 7 4];
colorsArea = cbrewer('qual', 'Set2', 8);
colorsArea = colorsArea(order,:);
colorsArea_sub = cbrewer('qual', 'Pastel2', 8);
colorsArea_sub = colorsArea_sub(order,:);
figure;
for p = 1 : length(test)
    clear bar_both bar_X bar_M
    eval(['bar_both = [sigSat_' test{p} '(:,1,1) sigSat_' test{p} '(:,2,1) ];'])
    eval(['bar_M = [sigSat_' test{p} '_M(:,1,1) sigSat_' test{p} '_M(:,2,1) ];'])
    eval(['bar_X = [sigSat_' test{p} '_X(:,1,1) sigSat_' test{p} '_X(:,2,1) ];'])
    subplot(1,length(test),p)
    for ar = 1 : length(area2test)
        bar(ar, (bar_both(ar,1)./bar_both(ar,2))*100,'FaceColor',colorsArea(ar,:));hold on
    end
    plot((bar_M(:,1)./bar_M(:,2))*100,'o','MarkerSize',10,'MarkerFaceColor',[.8 .8 .8],'MarkerEdgeColor','k')
    plot((bar_X(:,1)./bar_X(:,2))*100,'v','MarkerSize',10,'MarkerFaceColor',[.65 .65 .65],'MarkerEdgeColor','k')
    set(gca,'Xtick',1:length(area2test),'XtickLabel',area2test,'XtickLabelRotation',30,'FontSize',16)
    ylim([0 80])
    xlim([0 length(area2test)+1])
    ylabel('Percent neurons w/ sig Pref modulation')
    title([test{p} ' units'])
end

%- binomial model for perc Sig Neurons
x1= [];x2=[];x3 = [];
for mm = 1 : 2
    if mm ==1
        nN = sigSat_proba_M(:,:,1);
        % nN = sigSat_juice_M(:,:,1);
         nN = sigSat_smthg_M(:,:,1);
       %  nN = sigSat_all_M(:,:,1);
    else
        nN = sigSat_proba_X(:,:,1);
        % nN = sigSat_juice_X;
         nN = sigSat_smthg_X(:,:,1);
       %  nN = sigSat_all_X(:,:,1);
    end
    for i  = 1 : length(nN(:,1))
        x1 = [x1 ; repmat(area2test(i),nN(i,2),1)];
        x2 = [x2 ; repmat(true,nN(i,1),1) ; repmat(false,nN(i,2)-nN(i,1),1)];
    end
    x3 = [x3 ; repmat(['mk' num2str(mm)],sum(nN(:,2)),1)];
end

modeldata = table(logical(x2),x1,x3,'VariableNames',{'sig' 'area' 'mk'})

lme_percSig = fitglme(modeldata,'sig ~ 1 + area + (1|mk) ','Distribution','binomial','Link','logit')

anova(lme_percSig)
[~,wald_percNeurons,~,pval_adj_percNeurons] = area_posthoc(lme_percSig,area2test,'y');

test = {'proba' 'juice' 'smthg'};
clear chi2 pval
for t = 1 : length(test)
    for ar = 1 : length(area2test)
       eval(['[~,chi2(ar,t),pval(ar,t)] =  chi2_fms(sigSat_'  test{t} '(ar,1,1),sigSat_' test{t} '(ar,2,1),sigSat_nthg(ar,1,1),sigSat_nthg(ar,2,1));'])
    end
end
test = {'proba'};
clear chi2 pval
for t = 1 : length(test)
    for ar = 1 : length(area2test)
       eval(['[~,chi2(ar,t),pval(ar,t)] =  chi2_fms(sigSat_'  test{t} '(ar,1,1),sigSat_' test{t} '(ar,2,1),sigSat_juice(ar,1,1),sigSat_juice(ar,2,1));'])
    end
end

%- PLOT FOR WHETHER neurons with SAT more likely to decrease or increase their FR?
figure;
    for ar = 1 : length(area2test)
        bar(ar, (prop_SatIncr(ar,1)./prop_SatIncr(ar,2))*100,'FaceColor',colorsArea(ar,:));hold on
    end
    plot((prop_SatIncr_M(:,1)./prop_SatIncr_M(:,2))*100,'o','MarkerSize',10,'MarkerFaceColor',[.8 .8 .8],'MarkerEdgeColor','k')
    plot((prop_SatIncr_X(:,1)./prop_SatIncr_X(:,2))*100,'v','MarkerSize',10,'MarkerFaceColor',[.65 .65 .65],'MarkerEdgeColor','k')
    set(gca,'Xtick',1:length(area2test),'XtickLabel',area2test,'XtickLabelRotation',30,'FontSize',16)
    ylim([0 80])
    xlim([0 length(area2test)+1])
    ylabel('Percent of SAT neurons with POSITIVE BETAS')
    title([test{p} ' units'])












%- what's the most common motif
%- check against perm.
test = {'SAT' 'SATxPB' 'SATxJU' 'SATxPBxJU' };
clear pval_motif
for t = 1 : length(test)
    sigSat_true = satUnits(:,t,1)./satUnits(:,end,1);
    sigSat_perm = squeeze(satUnits(:,t,2:end))./squeeze(satUnits(:,end,2:end));
    pval_motif(:,t) = sum(sigSat_true<=sigSat_perm,2)/size(sigSat_perm,2); %- one-sided test
    rand_motif(:,t) = mean(sigSat_perm,2);
end

xax_ticks = ((0:1/(length(area2test)-1):1)-.5)*.675;
figure;
b = bar((satUnits(:,1:end-1,1)./repmat(satUnits(:,end,1),1,4))');hold on
for ar = 1 : length(b)
    b(ar).FaceColor = colorsArea(ar,:)
end
for t = 1 : length(test)
    plot(t+xax_ticks,rand_motif(:,t),'Color',[.6 .6 .6])
end 
for t = 1 : length(test)
    bar_M = satUnits_M(:,t,1)./satUnits_M(:,end,1);
    bar_X = satUnits_X(:,t,1)./satUnits_X(:,end,1);
    plot(t+xax_ticks,bar_M,'o','MarkerSize',10,'MarkerFaceColor',[.8 .8 .8],'MarkerEdgeColor','k')
    plot(t+xax_ticks,bar_X,'v','MarkerSize',10,'MarkerFaceColor',[.65 .65 .65],'MarkerEdgeColor','k')
end
ylabel('Percent of SAT neurons')
set(gca,'XTick',1:length(test),'XTickLabel',test,'FontSize',16)


%% plot examples

% figure;
% next = 0;
% for s = 1 : length(res_anova)
%     for n = 1 : length(res_anova(s).Ap(1,:))
%         if sum(res_anova(s).Ap([3 5 6 7],n)<0.001)~=0
%             next =next+1;
%             if next>3*15
%                 next = 1;figure;
%             end
%             subplot(3,15,next)
%                   params = {'ChPb' 'Ju' 'Satiety' 'ChPb x Ju' 'ChPb x Sat' ' Ju x Sat' 'ChPb x Ju x Sat'}
%                    
%                     colo = cbrewer('qual','Paired',6); colo = colo([1 2 5 6],:);
%                     xx=0;
%                    [a,b,c]= grpstats(res_anova(s).FR(n,:),res_anova(s).factors,{'mean' 'sem' 'gname'});
%                    c = cellfun(@str2num,c);
%                     nb_ju = length(unique(c(:,2)));
%                     nb_sat = length(unique(c(:,3)));
%                     for a1 = 1 : nb_ju
%                         for a2 = 1 : nb_sat
%                            idx=c(:,2)==a1 & c(:,3)==a2 
%                            xx=xx+1;
%                             plot(c(idx,1),a(idx),'.-','LineWidth',a2,'Color',colo(xx,:),'MarkerSize',25); hold on
%                            ciplot(a(idx)-b(idx),a(idx)+b(idx),c(idx,1),colo(xx,:),0.2)
%                         end
%                         
%                     end
%                     title([res_anova(s).neurons_area(n) params(res_anova(s).Ap(:,n)<0.001)])
%                     
%                 end
% 
% 
%     end
% end

                
%% plot examples PSTH

%- for saving some examples..
pathFig = 'C:\Users\Fred\Dropbox\Rudebeck Lab\Posters and Talks\PSTH Pref Example Neurons\';

for s = 1 : length(res_anova)
    for n = 1 : length(res_anova(s).Ap(1,:))
        clearvars -except s n res_anova pathFig
        if sum(res_anova(s).Ap([3 5 6 7],n)<0.0001)~=0
   %     if res_anova(s).Ap(6,n)<0.01

            %- extract the name of the file where the timestamps are!
            if res_anova(s).neurons_info(n,1)<10; addzeros = '00';
            elseif res_anova(s).neurons_info(n,1)<100; addzeros = '0';
            else addzeros = '';
            end
            if res_anova(s).neurons_info(n,2)<10; addzeros2 = '0';
            else addzeros2 = '';
            end

            file_id = [res_anova(s).session '_Ch' addzeros num2str(res_anova(s).neurons_info(n,1)) ...
                '_Clus' addzeros2 num2str(res_anova(s).neurons_info(n,2)) '.mat'];

            %- find the file in the different folder
            if strcmp(res_anova(s).session(1),'M')
                path4file = utils_POTT_SPKfolder('MORBIER');
            elseif strcmp(res_anova(s).session(1),'X') & datenum(res_anova(s).session(2:7),'mmddyy')<datenum('090920','mmddyy')
                path4file = utils_POTT_SPKfolder('MIMIC1');
            else
                path4file = utils_POTT_SPKfolder('MIMIC2');
            end

            pref_ch = res_anova(s).pref_ch;

            %- load unit
            Unit = load([path4file file_id]);
            Unit.info.area{1} = res_anova(s).neurons_area{n}; %- update area label based on histo

            %- load the event file from that session
            evt = load([path4file res_anova(s).session '_behav.mat']);

            %- take only the event considered in the Pref analysis
            TrialType_INS = evt.TrialType(:,evt.TrialType(1,:)==2 & evt.TrialType(2,:)==0 );

            evt2take = evt.t_stim(1,evt.TrialType(1,:)==2 & evt.TrialType(2,:)==0 ) ;

            dumm = TrialType_INS(ismember(evt.TrialType_header,{ 'I_chosenproba' 'I_chosenjuice' }),:)';
            dumm = [(1:length(dumm))' dumm(:,1)+(100*dumm(:,2))];

            when = @(arg,cd) TrialType_INS(strcmp(evt.TrialType_header,arg),:)==cd ;
            diff_juice = (when('I_juiceL',1) & when('I_juiceR',2)) | (when('I_juiceL',2) & when('I_juiceR',1)); %- take only the different juice trials

            %- extract only the diff juice trials
            tt = find(diff_juice);

            keep = [130 150 170 190 230 250 270 290];
            %  remove = ~ismember(dumm(:,2),keep) | ~diff_juice';

            preTr = pref_ch{1,2}; %- out of the diff juice trials
            preTr4spk = tt(preTr); %- trial to take for neurons
            preTr_fact = dumm(preTr4spk,:);
            preRemove = ~ismember(preTr_fact(:,2),keep) ;
            preTr_fact(preRemove,:)=[];

            postTr = pref_ch{1,3}; %- out of the diff juice trials
            postTr4spk = tt(postTr); %- trial to take for neurons
            postTr_fact = dumm(postTr4spk,:);
            postRemove = ~ismember(postTr_fact(:,2),keep) ;
            postTr_fact(postRemove,:)=[];

            if pref_ch{1,4}(2)-pref_ch{1,4}(1)>0
                allfactors = [preTr_fact(:,1) (mod(preTr_fact(:,2),100)) abs(floor(preTr_fact(:,2)/100)-3) ones(size(preTr_fact(:,1))) ; ...
                    postTr_fact(:,1) (mod(postTr_fact(:,2),100)) abs(floor(postTr_fact(:,2)/100)-3) 2*ones(size(postTr_fact(:,1)))]; %- inverse juice (1 = valued / 2 = devalued)
            else
                allfactors = [preTr_fact(:,1) (mod(preTr_fact(:,2),100)) floor(preTr_fact(:,2)/100) ones(size(preTr_fact(:,1))) ; ...
                    postTr_fact(:,1) (mod(postTr_fact(:,2),100)) floor(postTr_fact(:,2)/100) 2*ones(size(postTr_fact(:,1)))]; %- no inverse juice (1 = valued / 2 = devalued)
            end

            %- create conditions and threshold for subplot
            cond=allfactors(:,2)+100*allfactors(:,3)+1000*allfactors(:,4);
            evt2take = evt2take(allfactors(:,1));
            grp=cond>2000;

            %- general params
            evtName = 'Stim';
            pre = 750;
            post = 2000;

            fighandle1 = rasterPSTH_2bins(Unit,evt2take,evtName,cond',2000,pre,post,pathFig);
         %   fighandle1 = rasterPSTH_2bins(Unit,evt2take,evtName,cond',2000,pre,post);
         %   pause
            close(fighandle1)
        end
    end
end

                
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
    













