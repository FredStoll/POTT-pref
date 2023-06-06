%% Correlation of FR and popoulation activity with preference change

clear

monkey = 'X'

path2go = '/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-BehavPrefChange/data/neurons/subset-final/'; %- path where SPKpool files are!
path2go2 = '/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-BehavPrefChange/data/';
if strcmp(monkey,'X')
    load([path2go2 'Mimic_behav_bins_perms.mat'])
    list = dir([path2go 'X*a_SPKpool.mat']);
else
    load([path2go2 'Morbier_behav_bins_perms.mat'])
    list = dir([path2go 'M*a_SPKpool.mat']);
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
area2test = {'vlPFC' 'OFC' 'IFG' 'LAI' 'AMG'};

minNeurons = 5;
ref_length = 10;
normal = true;

meannorm  = @(data) (data-mean(mean(data)))/(max(max(data))-min(min(data)));
minnorm  = @(data) (data-min(min(data)))/(max(max(data))-min(min(data)));

pref_sig = NaN(length(list),4);
all_data=cell(2,length(list));
all_data_pop=cell(length(area2test),length(list));

allP = [];
allEst = [];
allID = [];
allCorr = [];
allCorr_full =[];
allP_pop = [];
allEst_pop = [];
allCorr_pop = [];
allCorr_full_pop =[];

for s =  1 : length(list)
    %for s =  [75 76 81 85 89]
    disp(s)
    SPK = load([path2go list(s).name]);
    idxSess = ismember(days,list(s).name(1:8));

    if ~isempty(SPK.neurons_area) %- skip sessions where no recorded neurons where in area2test!
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

        bins = param(idxSess).binSize;
        bins_start = round(1:(height(ALL(idxSess).T)-bins)/50:height(ALL(idxSess).T)-bins);

        %- find trial ID included for each 50 bins in the behavioral model
        clear trID
        dumm = find(ALL(idxSess).trial2take);
        for b = 1 : length(bins_start)
            trID(b,:) = dumm(bins_start(b):bins_start(b)+bins);
        end

        % find trial ID for the neurons
        trID_neurons = find(ALL(idxSess).TrialType(1,:)==2 & ALL(idxSess).TrialType(2,:)==0);

        [ALL(idxSess).Z_trend_ft,~]=Mann_Kendall(ALL(idxSess).ft_bins,0.01);
        [ALL(idxSess).Z_trend_rt,~]=Mann_Kendall(ALL(idxSess).rt_bins,0.01);
        [ALL(idxSess).Z_trend_pref,~]=Mann_Kendall(ALL(idxSess).pref_bins,0.01);

        if length(trID_neurons)~=size(SPK.TrialType_INS,2)
            warning('Problem in # of trials......')
            x = x + 1;
            pb_files(x)=s;

        else

            %- area for each neurons (1 to N)
            spk_area = zeros(size(SPK.neurons_area));
            for ar = 1 : length(area2test)
                eval(['dumm = find(ismember(SPK.neurons_area,area_list.' area2test{ar} '));'])
                spk_area(dumm) = ar;
            end

            unit_id = unique(SPK.Tr_Clust_INS(:,2));

            clear unit_fr
            for n = 1 : length(unit_id)
                for b = 1 : size(trID,1)
                    tr2use = find(ismember(trID_neurons,trID(b,:))==1); %- use the average firing rate across trials in each bin
                    unit_fr(n,b,:) = mean(mean(SPK.SPK_INS(ismember(SPK.Tr_Clust_INS(:,1),tr2use) & SPK.Tr_Clust_INS(:,2)==unit_id(n),:)));
                end
                %- normalization
                if normal && sum(sum(unit_fr(n,:)))~=0
                    unit_fr(n,:) = minnorm(unit_fr(n,:) );
                end
            end

            %  XX = squeeze(mean(unit_fr(n,:,:),3))';  %- extract new data (all trials and all times)
            XX = unit_fr';
            XX_ref = XX-repmat(mean(XX(1:ref_length,:)),size(XX,1),1);

            for pp = 1 %: length(ALLperm(1,:))+1
                if pp == 1
                    pref = ALL(idxSess).pref_bins-mean(ALL(idxSess).pref_bins(1:ref_length));
                    ft = ALL(idxSess).ft_bins-mean(ALL(idxSess).ft_bins(1:ref_length));
                    rt = ALL(idxSess).rt_bins-mean(ALL(idxSess).rt_bins(1:ref_length));
                else
                    pref = ALLperm(idxSess,pp-1).pref_bins-mean(ALLperm(idxSess,pp-1).pref_bins(1:ref_length));
                    ft = ALLperm(idxSess,pp-1).ft_bins-mean(ALLperm(idxSess,pp-1).ft_bins(1:ref_length));
                    rt = ALLperm(idxSess,pp-1).rt_bins-mean(ALLperm(idxSess,pp-1).rt_bins(1:ref_length));
                end

                %- FOR SINGLE NEURON CORRELATIONS
                all_res = []; all_r_temp = [];
                for n = 1 : length(unit_id)

                    mdl = fitglm([ft' , rt'],XX_ref(:,n),'VarNames',{'IT' 'RT' 'spk'});
                    allP = [allP ; mdl.Coefficients.pValue(2:end)' , spk_area(n) , s , pp-1];
                    allEst = [allEst ; mdl.Coefficients.Estimate(2:end)' , spk_area(n) , s , pp-1];

                    [r,p]=corrcoef(mdl.Residuals.Raw,pref);
                    allCorr = [allCorr ; r(2,1) , p(2,1) , spk_area(n) , s , pp-1];

                    %- for plotting purposes
                    all_r_temp = [all_r_temp ; r(2,1) p(2,1)];
                    all_res = [all_res , mdl.Residuals.Raw];

                    [r,p]=corrcoef(XX_ref(:,n),pref);corr_RP(:,1) = [r(2,1) ; p(2,1) ; ];
                    allCorr_full = [allCorr_full ; r(2,1) , p(2,1) , spk_area(n) , s , pp-1];

                end
                if pp == 1
                    all_data{1,s}=all_res;
                    all_data{2,s}=XX_ref;
                    allID=[allID ; SPK.neurons];
                end

                %- for plotting example
%                 nSpk = randperm(length(spk_area),10);
%                 data_plot = unit_fr(nSpk,:) + repmat([0:.3:.3*(length(unit_fr(nSpk,1))-1)]',1,50);
%                 datares_plot = all_res(:,nSpk)' + repmat([0:.3:.3*(length(all_res(1,nSpk)')-1)]',1,50);
%                 spk_area_sub = spk_area(nSpk);
% 
%                 order = [3 1 2 5 7 4];
%                 colorsArea = cbrewer('qual', 'Set2', 8);
%                 colorsArea = colorsArea(order,:);
%                 colorsArea_sub = cbrewer('qual', 'Pastel2', 8);
%                 colorsArea_sub = colorsArea_sub(order,:);
% 
%                 figure;
%                 subplot(1,2,1)
%                 for i = 1 : length(data_plot(:,1))
%                     plot(data_plot(i,:),'Color',colorsArea(spk_area_sub(i),:));hold on
%                 end
%                 subplot(1,2,2)
%                 for i = 1 : length(datares_plot(:,1))
%                     plot(datares_plot(i,:),'Color',colorsArea(spk_area_sub(i),:));hold on
%                 end

                %- FOR POPULATION DISTANCE
                for ar = 1 : length(area2test)
                    nNeurons(ar) = sum(spk_area==ar);

                    if nNeurons(ar)>=minNeurons

                      %  XX = unit_fr(spk_area==ar,:)';
                      %  XX_ref = XX-repmat(mean(XX(1:ref_length,:)),size(XX,1),1);
                        XX_ref = all_res(:,spk_area==ar);
                        [coeff,score,~,~,explained,~] = pca(XX_ref,'NumComponents',3);

                        spk_ref = mean(score(1:ref_length,1:3));
                        spk_distance = [];
                        for i = 1 : size(score,1)
                            spk_distance(i) = sqrt(sum((score(i,1:3)-spk_ref).^2));
                        end

                        %- for example
%                         figure;plot3(score(:,1),score(:,2),score(:,3),'Color',colorsArea(ar,:));hold on; plot3(spk_ref(1),spk_ref(2),spk_ref(3),'.','Color',colorsArea(ar,:),'MarkerSize',30)
%                         figure;plot(spk_distance,'Color',colorsArea(ar,:))
                     
                        [r,p]=corrcoef(spk_distance,pref);
                        allCorr_pop = [allCorr_pop ; r(2,1) , p(2,1) , ar , s , pp-1];

                        if pp == 1
                            all_data_pop{ar,s}=[score , spk_distance' , NaN(size(spk_distance')) , pref' , ft' , rt'];
                        end


                        XX = unit_fr(spk_area==ar,:)';
                        XX_ref = XX-repmat(mean(XX(1:ref_length,:)),size(XX,1),1);
                        [coeff,score,~,~,explained,~] = pca(XX_ref,'NumComponents',3);

                        spk_ref = mean(score(1:ref_length,1:3));
                        spk_distance = [];
                        for i = 1 : size(score,1)
                            spk_distance(i) = sqrt(sum((score(i,1:3)-spk_ref).^2));
                        end

                        [r,p]=corrcoef(spk_distance,pref);
                        allCorr_full_pop = [allCorr_full_pop ; r(2,1) , p(2,1) , ar , s , pp-1];

                       
                    end
                end


            end
            pref_sig(s,:) = [ALL(idxSess).p_trend ALL(idxSess).Z_trend_pref ALL(idxSess).Z_trend_ft ALL(idxSess).Z_trend_rt ];

        end
    end
end

save([path2go2 monkey '_Corr_pref_STIM_final.mat'],'pref_sig','all_data','allCorr','allCorr_full','allEst','allP','allID',...
    'all_data_pop','allCorr_pop','allCorr_full_pop','allEst_pop','allP_pop','area2test')            
%save([path2go2 monkey '_Corr_pref_REF.mat'],'pref_sig','all_unit_data','allCorr','allEst','allP','area2test')            

%% 

clear
path2go2 = '/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-BehavPrefChange/data/';
load([path2go2 'M_Corr_pref_STIM_final.mat']);            
X = load([path2go2 'X_Corr_pref_STIM_final.mat']);      

Stability = load('C:\Users\Fred\Dropbox\Rudebeck Lab\ANA-POTT-BehavPrefChange\data\POTT_waveforms_ratio.mat');
rej = false; %- reject if waveform not fully stable

%- COLOR ASSIGNMENT

order = [3 1 2 5 4 7];
colorsArea = cbrewer('qual', 'Set2', 8);
colorsArea = colorsArea(order,:);
colorsArea_sub = cbrewer('qual', 'Pastel2', 8);
colorsArea_sub = colorsArea_sub(order,:);


%- combine monkeys
mk = [zeros(size(pref_sig,1),1) ;ones(size(X.pref_sig,1),1)  ];
pref_sig = [pref_sig , (1:length(pref_sig))' ; X.pref_sig 1000+(1:length(X.pref_sig))'];
allCorr = [allCorr ; [X.allCorr(:,1:3) , X.allCorr(:,4)+1000 , X.allCorr(:,5)]];
allCorr_full = [allCorr_full ; [X.allCorr_full(:,1:3) , X.allCorr_full(:,4)+1000 , X.allCorr_full(:,5)]];
allCorr_pop = [allCorr_pop ; [X.allCorr_pop(:,1:3) , X.allCorr_pop(:,4)+1000 , X.allCorr_pop(:,5)]];
allCorr_full_pop = [allCorr_full_pop ; [X.allCorr_full_pop(:,1:3) , X.allCorr_full_pop(:,4)+1000 , X.allCorr_full_pop(:,5)]];
allID = [allID ; X.allID ];

[h_sig, crit_p, adj_p]=fdr_bh(pref_sig(:,1),0.05,'pdep','yes');
pref_sig(:,1) = adj_p;

%allCorr = allCorr_full;
%allCorr_pop = allCorr_full_pop;

%- find the sig sessions for pref
for i = 1 : length(allCorr)
    xx = find(pref_sig(:,5)==allCorr(i,4));
    allCorr(i,5:6)= [pref_sig(xx,1) mk(xx)];
end
for i = 1 : length(allCorr_pop)
    xx = find(pref_sig(:,5)==allCorr_pop(i,4));
    allCorr_pop(i,5:6)= [pref_sig(xx,1) mk(xx)];
end

%- match isolation measures
stab_measures=[];
for n = 1 : length(allID) 
    loc2match =  find(ismember(Stability.spk_names,{allID(n,:)}));
    stab_measures(n,:)=Stability.avg_ratio_sig(loc2match,:) ;

%    [Z_trend_wf,P_wf]=Mann_Kendall(Stability.avg_ratio_norm(loc2match,:),0.01);
%    stab_measures_pval(n,1) = P_wf;
end

all_data = [all_data , X.all_data];
all_data_pop = [all_data_pop , X.all_data_pop];


if rej
    reject =  (sum(stab_measures')~=0)';
   % reject =  stab_measures_pval<0.01;
    allCorr(isnan(allCorr(:,2)) | reject,:)=[];
    allCorr_full(isnan(allCorr_full(:,2)) | reject,:)=[];
else
    allCorr(isnan(allCorr(:,2)) ,:)=[];
    allCorr_full(isnan(allCorr_full(:,2)) ,:)=[];
end

allCorr_pop(isnan(allCorr_pop(:,2)),:)=[];
allCorr_full_pop(isnan(allCorr_full_pop(:,2)),:)=[];

%- ignore that cos no permutation for now
% allCorr_all = allCorr;
% allCorr = allCorr(allCorr(:,end)==0,1:4);




% allR =[];
% figure;
% for ar = 1 : length(area2test)
%     [rr,pp]=corrcoef(allCorr(allCorr(:,3)==ar &  allCorr(:,5)<0.05,2),allCorr_full(allCorr(:,3)==ar &  allCorr(:,5)<0.05,2))
%     allR(ar,:) = [rr(2,1) pp(2,1)];
%     subplot(1,length(area2test),ar);
%     plot(allCorr(allCorr(:,3)==ar &  allCorr(:,5)<0.05,1),allCorr_full(allCorr(:,3)==ar &  allCorr(:,5)<0.05,1),'o')
% end
%- reject AMG, too little variability in sessions
allCorr(allCorr(:,3)==find(ismember(area2test,'AMG')),:)=[];
allCorr_pop(allCorr_pop(:,3)==find(ismember(area2test,'AMG')),:)=[];


sigPref = allCorr(:,5)<0.05;
sigPref_pop = allCorr_pop(:,5)<0.05;

measures = {'NEURONS' 'SESSIONS'}
for m =  1 : 2

    if m==1
        modeldata = table(double(abs(allCorr(sigPref,1))),area2test(allCorr(sigPref,3))',categorical(allCorr(sigPref,6)),categorical(allCorr(sigPref,4)), 'VariableNames',{'perf' 'area' 'mk' 'sess'})
    else
        modeldata = table(double(abs(allCorr_pop(sigPref_pop,1))),area2test(allCorr_pop(sigPref_pop,3))',categorical(allCorr_pop(sigPref_pop,6)),categorical(allCorr_pop(sigPref_pop,4)), 'VariableNames',{'perf' 'area' 'mk' 'sess'})
    end 
    
    lme = fitglme(modeldata,'perf ~ 1 + area  + (1|mk) + (1|sess)')
    [pval,wald,thr_corr,pval_adj] = area_posthoc(lme,area2test(1:4),'y');
    anova(lme)
    pause;
    % grpstats(modeldata(:,1:3),{'area' 'mk'},"mean")
    
     
    
    [h, crit_p, adj_p]=fdr_bh(allCorr(:,2),0.05,'pdep','yes');
    [h_pop, crit_p, adj_p]=fdr_bh(allCorr_pop(:,2),0.05,'pdep','yes');
    
    if m == 1 
        sig = false(length(allCorr),1);
        sig(h,1)=true;
        modeldata = table(sig(sigPref),area2test(allCorr(sigPref,3))',categorical(allCorr(sigPref,6)),categorical(allCorr(sigPref,4)), 'VariableNames',{'sig' 'area' 'mk' 'sess'})
    else
        sig = false(length(allCorr_pop),1);
        sig(h_pop,1)=true;
        modeldata = table(sig(sigPref_pop),area2test(allCorr_pop(sigPref_pop,3))',categorical(allCorr_pop(sigPref_pop,6)),categorical(allCorr_pop(sigPref_pop,4)), 'VariableNames',{'sig' 'area' 'mk' 'sess'})
    end
           
           lme = fitglme(modeldata,'sig ~ 1 + area  + (1|mk) + (1|sess)','Distribution','Binomial')

            [pval,wald,thr_corr,pval_adj] = area_posthoc(lme,area2test(1:4),'y');
           anova(lme)
           disp([sum(modeldata.sig) length(modeldata.sig)])
           pause;
    
    
           for ar = 1 : length(area2test)
               nbSig(ar,:) = [sum(modeldata.sig(ismember(modeldata.area,area2test{ar})) ) length(modeldata.sig(ismember(modeldata.area,area2test{ar})) )];
               dumm = modeldata.sig(ismember(modeldata.area,area2test{ar}) & modeldata.mk==categorical(0));
               nbSig_M(ar,:) = [sum(dumm ) length(dumm)];
               dumm = modeldata.sig(ismember(modeldata.area,area2test{ar}) & modeldata.mk==categorical(1));
               nbSig_X(ar,:) = [sum(dumm ) length(dumm)];
           end
    
    
    figure;
    
    subplot(1,2,1)
        for ar = 1 : length(area2test)
            bar(ar, (nbSig(ar,1)./nbSig(ar,2))*100,'FaceColor',colorsArea(ar,:));hold on
        end
        plot((nbSig_M(:,1)./nbSig_M(:,2))*100,'o','MarkerSize',10,'MarkerFaceColor',[.8 .8 .8],'MarkerEdgeColor','k')
        plot((nbSig_X(:,1)./nbSig_X(:,2))*100,'v','MarkerSize',10,'MarkerFaceColor',[.65 .65 .65],'MarkerEdgeColor','k')
        set(gca,'Xtick',1:length(area2test),'XtickLabel',area2test,'XtickLabelRotation',30,'FontSize',16)
        ylim([0 100])
        xlim([0 length(area2test)+1])
        ylabel(['Percent ' measures{m} ' w/ sig Pref correlation'])
    
    subplot(1,2,2)
    for ar = 1 : length(area2test)-1
    
        yl=0+ar;
        wdth = .35;
    if m == 1 
        X = abs(allCorr(allCorr(:,3)==ar & sigPref ,1));
    else
        X = abs(allCorr_pop(allCorr_pop(:,3)==ar & sigPref_pop ,1));
    end
        boxplot_ind(X,yl,wdth,[colorsArea_sub(ar,:) ;colorsArea(ar,:)]);hold on
        set(gca,'view',[90 -90],'color','none','FontSize',16);
        
    end
    xlim([0 1]);ylim([0 7])
    xlabel('Correlation values')
    set(gca,'Ytick',1:length(area2test),'YtickLabel',area2test,'YtickLabelRotation',30,'FontSize',16)

end


%- Pairwise correlation difference

diff_corr=cell(6,6)
sess = unique(allCorr_pop(:,4))
for s = 1 : length(sess)
    mat = allCorr_pop(allCorr_pop(:,4)==sess(s),:);
    for ar1 = 1: length(area2test)
        for ar2 = 1: length(area2test)
            if ar1 ~= ar2 & ar1 < ar2 & sum(ismember(mat(:,3),[ar1 ar2]))==2
               diff_corr{ar1,ar2} = [diff_corr{ar1,ar2} ; abs(mat(mat(:,3)==ar1,1)) abs(mat(mat(:,3)==ar2,1))];
            end
        end
    end
end

diff_CC = [];
name_CC = {};
x=0;
    for ar1 = 1: length(area2test)
        for ar2 = 1: length(area2test)
            if ar1 ~= ar2 & ar1 < ar2 & ~isempty(diff_corr{ar1,ar2})
                x = x + 1;
                diff_CC(1,x) = mean(diff_corr{ar1, ar2}(:,2)-diff_corr{ar1, ar2}(:,1));
                diff_CC(2,x) = std(diff_corr{ar1, ar2}(:,2)-diff_corr{ar1, ar2}(:,1))/sqrt(length(diff_corr{ar1, ar2}));
                name_CC{x} = [area2test{ar2} '-' area2test{ar1}];
                dumm = [diff_corr{ar1, ar2}(:,2) , ones(length(diff_corr{ar1, ar2}(:,2)),1)  ; ...
                        diff_corr{ar1, ar2}(:,1) , 2*ones(length(diff_corr{ar1, ar2}(:,2)),1)]
                [p_val(x),b,c]=kruskalwallis(dumm(:,1),dumm(:,2),'off');

            end
        end
    end

figure;
line([0 0],[0 length(name_CC)+1],'Color','k');hold on;
plot(diff_CC(1,:),1:length(name_CC),'.','MarkerSize',20);
hold on;
for i = 1 : length(name_CC)
    line([diff_CC(1,i)-diff_CC(2,i) diff_CC(1,i)+diff_CC(2,i)],[i i],'Color','k')
end
set(gca,'YTick',[1:length(name_CC)],'YTickLabel',name_CC)
xlim([-.15 .15]);box on
xlabel('Rho difference (ar1-ar2)')


%- example figure

for sess = 182 % 1 : length(all_data_pop(1,:));
    pop =[];
    ars=[];
    pref=[];
    for ar = 1 : length(area2test)-1
        if ~isempty(all_data_pop{ar,sess})
            pref = all_data_pop{ar,sess}(:,6);
            pop = [pop all_data_pop{ar,sess}(:,4)];
            ars=[ars ar];
        end
    end
    if ~isempty(pref) & pref_sig(sess,1)<0.01
        f = figure('Position',[650   258   682   692]);

        subplot(4,6,[1:5]);
        plot( pref ,'color','k','LineWidth',2);xlim([.5 50.5]);%ylim([0 1])
        ylabel('Preference')
        title(['sess=' num2str(sess) ' - p=' num2str(pref_sig(sess,1))])
        set(gca,"FontSize",16)

        subplot(4,6,[7:11 13:17]);
        corrvalues = allCorr(allCorr(:,4) == pref_sig(sess,end),[1 3]);
        [a,b] = sortrows(corrvalues(:,1));
        imagesc(all_data{1,sess}(:,b)'); axis xy ;xlim([.5 50.5])
        ylabel('Neurons')
        title('GLM residuals - Neurons')
        set(gca,"FontSize",16)

        subplot(4,6,[12 18]);
        for ar = 1 : length(area2test)
            if sum(corrvalues(:,2)==ar)~=0
                dumm = corrvalues ;
                dumm(b(corrvalues(b,2)~=ar),1)=0;
                bb = barh( dumm(b,1)  ) ;
                bb.FaceColor = colorsArea(ar,:);
                hold on
            end
        end
        line([-1 1],[.5 .5],'Color','k')
        %  barh(all_r_temp(b,1)) ; hold on
        %  barh(find(sig(b)==1),all_r_temp(b(sig(b)),1)) ; hold on
        ylim([.5 length(b)+.5]);
        xlim([-1 1]); axis off
        xlabel('R')
        set(gca,"FontSize",16)

        subplot(4,6,19:23);
        for ar = 1 : length(ars)
            plot(pop(:,ar),'Color',colorsArea(ars(ar),:),'LineWidth',2);  hold on
            xlim([.5 50.5])
        end
        xlabel('Time bins')
        title('GLM residuals - Populations')
        set(gca,"FontSize",16)

     %   pause;
     %   close(f)
    end
end



