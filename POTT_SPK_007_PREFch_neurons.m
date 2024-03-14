%% POTT Pref - SPK ANOVA for preference modulation - Figure 7
%- Require the time windows for PRE/POST pref change extracted in
%- POTT_BHV_001_GLM (Matrix Pref_bins.mat)
%-
%- Author: Fred M. Stoll, Icahn School of Medicine at Mount Sinai, NY
%- Date: 2022.12
%- Related to: Stoll & Rudebeck, Neuron, 2024

%% Run the ANOVA on neurons' firing rates..
%- do that for each monkey (X and M)

clear

monkey = 'X'; %- X / M
norm_me = @(data) -1+((data-min(data))*2)/(max(data)-min(data)) ;
norm_pb = @(data) -1+((data-0.1)*2)/(0.9-0.1) ;

path2go = 'C:\\Users\fred\Dropbox\Rudebeck Lab\ANA-POTT-BehavPrefChange\data\neurons\subset-final\';
path2go2 = 'C:\\Users\fred\Dropbox\Rudebeck Lab\ANA-POTT-BehavPrefChange\data\';

load([path2go2 'Pref_bins.mat'])

if strcmp(monkey,'X')
    load([path2go2 'Mimic_behav_bins.mat'])
    list = dir([path2go 'X*a_SPKpool.mat']);
    pref_ch = pref_ch_mk{2};
else
    load([path2go2 'Morbier_behav_bins.mat'])
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
bins4decoding=[2]; %- perform decoding on subset on bins (stim period here)
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
        time2cons = ismember(bins,bins4decoding) & time>=binwindow(1) & time<=binwindow(2);
        SPK.SPK_INS = SPK.SPK_INS(:,time2cons);
        time(~time2cons) = [];
              
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
                        
            dumm = SPK.TrialType_INS(ismember(SPK.TrialType_header,{'nTr' 'I_chosenproba' 'I_chosenjuice' }),:)';
            dumm = [dumm(:,1) dumm(:,2)+(100*dumm(:,3))];
            
            when = @(arg,cd) SPK.TrialType_INS(strcmp(SPK.TrialType_header,arg),:)==cd ;
            diff_juice = (when('I_juiceL',1) & when('I_juiceR',2)) | (when('I_juiceL',2) & when('I_juiceR',1)); %- take only the different juice trials
            
            %- extract only the diff juice trials
            tt = find(diff_juice);
            
            keep = [130 150 170 190 230 250 270 290];
            
            preTr = pref_ch{takebhv,2}; %- out of the diff juice trials
            preTr4spk = tt(preTr); %- trial to take for neurons PRE
            preTr_fact = dumm(preTr4spk,:);
            preRemove = ~ismember(preTr_fact(:,2),keep) ;
            preTr_fact(preRemove,:)=[];
          
            postTr = pref_ch{takebhv,3}; %- out of the diff juice trials
            postTr4spk = tt(postTr); %- trial to take for neurons POST
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

save([path2go 'res_ANOVAprefbin_' monkey '_final.mat'],'res_anova','list','-v7.3')

%% Post hoc

clear
path2go = '/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-BehavPrefChange/data/neurons/subset-final/'; %- path where SPKpool files are!

Stability = load('C:\Users\Fred\Dropbox\Rudebeck Lab\ANA-POTT-BehavPrefChange\data\POTT_waveforms_ratio.mat')
rej = false; %- reject if waveform not fully stable

X = load([path2go 'res_ANOVAprefbin_X_final.mat'])
M = load([path2go 'res_ANOVAprefbin_M_final.mat'])

res_anova = [X.res_anova M.res_anova];
area2test = {'vlPFC' 'OFC' 'IFG' 'LAI'};
area_list = utils_POTT_areas;

thr_sig = 0.01;

%- to show the stat on figures
statout = @(out,p2take_name) ['F(' num2str(out.DF1(strcmp(out.Term,p2take_name))) ',' num2str(out.DF2(strcmp(out.Term,p2take_name))) ...
    ')=' num2str(round(out.FStat(strcmp(out.Term,p2take_name)),3)) ...
    ', p(' p2take_name ')=' num2str(round(out.pValue(strcmp(out.Term,p2take_name)),3))];

%- combine all sessions
all_diff = [res_anova(:).Ap];
all_diff_perm = [res_anova(:).Ap_perm];
all_units = cat(1,res_anova(:).neurons_area);
all_beta = [res_anova(:).Betas];
all_omega = [res_anova(:).Omega2];
all_omega_perm = [res_anova(:).Omega2_perm];

all_mk = [];
for s = 1 : length(res_anova)
    all_mk = [all_mk ; repmat({res_anova(s).session(1)},size(res_anova(s).neurons_area,1),1)];
end

%- remove neurons with ANOVAs not converging
keep = (sum(isnan(all_diff))==0)' & ~ismember(all_units,'AMG');

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

%- session ID (random order)
all_sess = cat(1,res_anova(:).session);
all_sess = all_sess(randperm(length(all_sess)),:);
for n = 1 : length(all_names) 
    all_sess_id(n,1)=find(ismember(all_sess,{all_names(n,1:8)}));
end

%- match isolation measures
stab_measures=[];
for n = 1 : length(all_names) 
    loc2match =  find(ismember(Stability.spk_names,{[all_names(n,:) '.mat']}));
    stab_measures(n,:)=Stability.avg_ratio_sig(loc2match,:) ;
end

%- if removing some units depending on isolation measures...
if rej
    keep = (sum(isnan(all_diff))==0)' & ~ismember(all_units,'AMG') & (sum(stab_measures')==0)' ;
end

probaUnits = sum(all_diff([1 4 5 7],keep)<=thr_sig)~=0 ;
juiceUnits = sum(all_diff([2 4 6 7],keep)<=thr_sig)~=0 ;
bothUnits = juiceUnits & probaUnits ;
probaonlyUnits = ~juiceUnits & probaUnits ;
juiceonlyUnits = juiceUnits & ~probaUnits ;
eitherUnits = juiceUnits | probaUnits ;
nthgUnits = ~juiceUnits & ~probaUnits ;
allUnits = true(size(juiceUnits)) ;

sig_units = sum(all_diff([3 5 6 7],keep)<=thr_sig)~=0;
[sum(sig_units) length(sig_units)]
sig_units_perm = squeeze(sum(all_diff_perm([3 5 6 7],keep,:)<=thr_sig))~=0;
[min(sum(sig_units_perm)) mean(sum(sig_units_perm)) max(sum(sig_units_perm))]
 
all_area = cell(size(all_units));
for ar = 1 : length(area2test)
    eval(['dumm = ismember(all_units,area_list.' area2test{ar} ') ;'])
    all_area(dumm) = area2test(ar) ;
end
all_area = all_area(keep);
all_mk = all_mk(keep);
all_beta = all_beta(:,keep);
all_sess_id = all_sess_id(keep);

%- count nb neurons
for ar = 1 : length(area2test)
    dumm = ismember(all_area, area2test{ar} ) ;
    nbunits_tot(ar,:) = [sum(dumm & ismember(all_mk,'M'))  sum(dumm & ismember(all_mk,'X')) ];
end

modeldata = table(sig_units',probaUnits',juiceUnits',all_area,all_mk,all_sess_id, 'VariableNames',{'sig' 'proba' 'juice' 'area' 'mk' 'sess'})
figure;
prop_venn=[];
colors4pie = [253 192 134 ; 190 174 212; 127 201 127 ; 153 153 153]/255;
for ar = 1 : length(area2test)
    sub = modeldata(ismember(modeldata.area,area2test{ar}) ,:);
    prop_venn(ar,:) = [sum(sub.sig==true & sub.proba==true & sub.juice==false) ...
         sum(sub.sig==true & sub.proba==true & sub.juice==true) ...
         sum(sub.sig==true & sub.proba==false & sub.juice==true) ...
         sum(sub.sig==true & sub.proba==false & sub.juice==false) ...
         ];

    prop = round(100*(prop_venn(ar,:)/sum(prop_venn(ar,:))))
    labels = {['PB - ' num2str(prop(1)) '%'],['PB+ID - ' num2str(prop(2)) '%'],['ID - ' num2str(prop(3)) '%'],['others - ' num2str(prop(4)) '%']};

    subplot(2,2,ar)
   % pat = pie(prop_venn(ar,:),labels);
    pat = pie(prop_venn(ar,:));
    x = 0;
    for i = 1 : 2 : length(pat)
        x = x + 1;
        pat(i).FaceColor = colors4pie(x,:);
    end
    title(area2test{ar})
end

modeldata = table(sig_units',probaUnits',juiceUnits',all_area,all_mk,all_sess_id, 'VariableNames',{'sig' 'proba' 'juice' 'area' 'mk' 'sess'});
modeldata = modeldata(modeldata.sig==true,:);
models_form = {'sig ~ 1 + area  + (1|mk) + (1|sess)'};

%- significance for proba only
cd = modeldata.proba==true & modeldata.juice==false;
modeldata_sub = table(cd, modeldata.area,modeldata.mk,modeldata.sess , 'VariableNames',{'sig' 'area' 'mk' 'sess'});
lme_pb = fitglme(modeldata_sub,models_form{1},'Distribution','Binomial');
[pval,wald,thr_corr,pval_adj] = area_posthoc(lme_pb,area2test,'n');
disp(['%%%%%%%%%% Proba only neurons with model = ' models_form{1} ' %%%%%%%%%%'])
disp(anova(lme_pb));disp(pval_adj);disp(wald);
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

%- significance for proba + flavor
cd = modeldata.proba==true & modeldata.juice==true;
modeldata_sub = table(cd, modeldata.area,modeldata.mk,modeldata.sess , 'VariableNames',{'sig' 'area' 'mk' 'sess'});
lme_pbfl = fitglme(modeldata_sub,models_form{1},'Distribution','Binomial'); model_final = models_form{1};
[pval,wald,thr_corr,pval_adj] = area_posthoc(lme_pbfl,area2test,'n');
disp(['%%%%%%%%%% Proba + Flavor neurons with model = ' model_final ' %%%%%%%%%%'])
disp(anova(lme_pbfl));disp(pval_adj);disp(wald);
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

%- significance for flavor
cd = modeldata.proba==false & modeldata.juice==true;
modeldata_sub = table(cd, modeldata.area,modeldata.mk,modeldata.sess , 'VariableNames',{'sig' 'area' 'mk' 'sess'});
lme_fl = fitglme(modeldata_sub,models_form{1},'Distribution','Binomial'); model_final = models_form{1};
[pval,wald,thr_corr,pval_adj] = area_posthoc(lme_fl,area2test,'n');
disp(['%%%%%%%%%% Proba + Flavor neurons with model = ' model_final ' %%%%%%%%%%'])
disp(anova(lme_fl));disp(pval_adj);disp(wald);
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

%- significance for others
cd = modeldata.proba==false & modeldata.juice==false;
modeldata_sub = table(cd, modeldata.area,modeldata.mk,modeldata.sess , 'VariableNames',{'sig' 'area' 'mk' 'sess'});
lme_ot = fitglme(modeldata_sub,models_form{1},'Distribution','Binomial'); model_final = models_form{1};
[pval,wald,thr_corr,pval_adj] = area_posthoc(lme_ot,area2test,'n');
disp(['%%%%%%%%%% Proba + Flavor neurons with model = ' model_final ' %%%%%%%%%%'])
disp(anova(lme_ot));disp(pval_adj);disp(wald);
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

%- show stat on Fig 7C
legend({statout(anova(lme_pb),'area') statout(anova(lme_pbfl),'area') statout(anova(lme_fl),'area') statout(anova(lme_ot),'area') })


%- for the different populations (stats on Fig 7D)
test = {'probaonly' 'both' 'juiceonly' 'nthg' 'all'};
for t = 1 : length(test)
    eval(['restr = ' test{t} 'Units;'])
        %- significance
        modeldata = table(sig_units(restr)',all_area(restr),all_mk(restr),all_sess_id(restr), 'VariableNames',{'sig' 'area' 'mk' 'sess'});
        lme = fitglme(modeldata,models_form{1}); 

    [pval,wald,thr_corr,pval_adj] = area_posthoc(lme,area2test,'n');
        disp(['%%%%%%%%%% % sig ' test{t} ' with model = ' model_final ' %%%%%%%%%%'])
        disp(anova(lme));disp(pval_adj);disp(wald);
        stats_fig{t} = statout(anova(lme),'area');
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

    for ar = 1 : length(area2test)
        nbunits_stab(ar,:) = [sum(modeldata.sig(ismember(modeldata.area,area2test{ar}))) length(modeldata.sig(ismember(modeldata.area,area2test{ar})))];
    end
    nbunits_stab
end

%- check permutations numbers
t=1;
for p = 1 : length(sig_units_perm(1,:))
    eval(['restr = ' test{t} 'Units;']);
    modeldata = table(sig_units_perm(restr,p),all_area(restr),all_mk(restr),all_sess_id(restr), 'VariableNames',{'sig' 'area' 'mk' 'sess'});

    for ar = 1 : length(area2test)
        nbSig(ar,:) = [sum(modeldata.sig(ismember(modeldata.area,area2test{ar})) ) length(modeldata.sig(ismember(modeldata.area,area2test{ar})) )];
        dumm = modeldata.sig(ismember(modeldata.area,area2test{ar}) & ismember(modeldata.mk,'M'));
        nbSig_M(ar,:) = [sum(dumm ) length(dumm)];
        dumm = modeldata.sig(ismember(modeldata.area,area2test{ar}) & ismember(modeldata.mk,'X'));
        nbSig_X(ar,:) = [sum(dumm ) length(dumm)];
    end
    nbSig_allperm(p,:)=sum(nbSig);
end
mean(nbSig_allperm(:,1)./nbSig_allperm(:,2))
min(nbSig_allperm(:,1)./nbSig_allperm(:,2))
max(nbSig_allperm(:,1)./nbSig_allperm(:,2))


%- COLOR ASSIGNMENT
order = [3 1 2 5 4 7];
colorsArea = cbrewer('qual', 'Set2', 8);
colorsArea = colorsArea(order,:);
colorsArea_sub = cbrewer('qual', 'Pastel2', 8);
colorsArea_sub = colorsArea_sub(order,:);

%- Figure 7D
figure;
for t = 1 : length(test)
    eval(['restr = ' test{t} 'Units;'])
    modeldata = table(sig_units(restr)',all_area(restr),all_mk(restr),all_sess_id(restr), 'VariableNames',{'sig' 'area' 'mk' 'sess'})

    for ar = 1 : length(area2test)
        nbSig(ar,:) = [sum(modeldata.sig(ismember(modeldata.area,area2test{ar})) ) length(modeldata.sig(ismember(modeldata.area,area2test{ar})) )];
        dumm = modeldata.sig(ismember(modeldata.area,area2test{ar}) & ismember(modeldata.mk,'M'));
        nbSig_M(ar,:) = [sum(dumm ) length(dumm)];
        dumm = modeldata.sig(ismember(modeldata.area,area2test{ar}) & ismember(modeldata.mk,'X'));
        nbSig_X(ar,:) = [sum(dumm) length(dumm)];
    end

    subplot(2,3,t)
    for ar = 1 : length(area2test)
        bar(ar, (nbSig(ar,1)./nbSig(ar,2))*100,'FaceColor',colorsArea(ar,:));hold on
    end
    plot((nbSig_M(:,1)./nbSig_M(:,2))*100,'o','MarkerSize',10,'MarkerFaceColor',[.8 .8 .8],'MarkerEdgeColor','k')
    plot((nbSig_X(:,1)./nbSig_X(:,2))*100,'v','MarkerSize',10,'MarkerFaceColor',[.65 .65 .65],'MarkerEdgeColor','k')
    set(gca,'Xtick',1:length(area2test),'XtickLabel',area2test,'XtickLabelRotation',30,'FontSize',16)
    ylim([0 100])
    xlim([0 length(area2test)+1])
    ylabel('Percent Pref modulated neurons')
    title({[test{t} ' units'],stats_fig{t}})
end

%% Check the sign of the Beta (decrease or increase in FR with PREF)
all_beta_sig = all_beta;
all_beta_sig(all_diff(:,keep)>thr_sig)=NaN;

for ar = 1 : length(area2test)
    dumm = sign(nanmean(all_beta_sig([3 5 6 7],ismember(all_area,area2test(ar)) & sig_units')));
    prop_SatIncr(ar,:) = [sum(dumm==1) length(dumm)];
    dumm = sign(nanmean(all_beta_sig([3 5 6 7],ismember(all_area,area2test(ar)) & sig_units' & ismember(all_mk,'M'))));
    prop_SatIncr_M(ar,:) = [sum(dumm==1) length(dumm)];
    dumm = sign(nanmean(all_beta_sig([3 5 6 7],ismember(all_area,area2test(ar)) & sig_units' & ismember(all_mk,'X'))));
    prop_SatIncr_X(ar,:) = [sum(dumm==1) length(dumm)];
end

figure;
for ar = 1 : length(area2test)
barh(ar, (prop_SatIncr(ar,1)./prop_SatIncr(ar,2))*100,'FaceColor',colorsArea(ar,:));hold on
end
plot((prop_SatIncr_M(:,1)./prop_SatIncr_M(:,2))*100,1:length(prop_SatIncr_M),'o','MarkerSize',10,'MarkerFaceColor',[.8 .8 .8],'MarkerEdgeColor','k')
plot((prop_SatIncr_X(:,1)./prop_SatIncr_X(:,2))*100,1:length(prop_SatIncr_X),'v','MarkerSize',10,'MarkerFaceColor',[.65 .65 .65],'MarkerEdgeColor','k')
set(gca,'Ytick',1:length(area2test),'YtickLabel',area2test,'FontSize',16)
xlim([0 100])
ylim([0 length(area2test)+1]);axis ij;
line([50 50],[0 length(area2test)+1],'Color','k')
xlabel('Percent of SAT neurons with POSITIVE BETAS')


prop_SatIncr_rand = [round(prop_SatIncr(:,2)/2) prop_SatIncr(:,2)];
for ar = 1 : length(area2test)
    [a,chi_bias(ar),p_bias(ar)]=chi2_fms(prop_SatIncr(ar,1),prop_SatIncr(ar,2),prop_SatIncr_rand(ar,1),prop_SatIncr_rand(ar,2));
end


%% what's the most common motif ? Figure S8
%- check against perm.
clear sat_units sat_units_M sat_units_X
all_diff_restr = all_diff(:,keep);
for ar = 1 : length(area2test)
    dumm_area = ismember(all_area,area2test(ar))
    sat_units(ar,:) = [sum(all_diff_restr([3],dumm_area)<=thr_sig  & sum(all_diff_restr([5 6 7],dumm_area)<=thr_sig)==0 ) ; ...
        sum(all_diff_restr([5],dumm_area)<=thr_sig  & sum(all_diff_restr([6 7],dumm_area)<=thr_sig)==0 ) ; ...
        sum(all_diff_restr([6],dumm_area)<=thr_sig  & sum(all_diff_restr([5 7],dumm_area)<=thr_sig)==0 ) ; ...
        sum(all_diff_restr([7],dumm_area)<=thr_sig) ; ...
        sum(sum(all_diff_restr([3 5 6 7],dumm_area)<=thr_sig)~=0)]';
    sat_units_M(ar,:) = [sum(all_diff_restr([3],dumm_area & ismember(all_mk,'M'))<=thr_sig  & sum(all_diff_restr([5 6 7],dumm_area & ismember(all_mk,'M'))<=thr_sig)==0 ) ; ...
        sum(all_diff_restr([5],dumm_area & ismember(all_mk,'M'))<=thr_sig  & sum(all_diff_restr([6 7],dumm_area & ismember(all_mk,'M'))<=thr_sig)==0 ) ; ...
        sum(all_diff_restr([6],dumm_area & ismember(all_mk,'M'))<=thr_sig  & sum(all_diff_restr([5 7],dumm_area & ismember(all_mk,'M'))<=thr_sig)==0 ) ; ...
        sum(all_diff_restr([7],dumm_area & ismember(all_mk,'M'))<=thr_sig) ; ...
        sum(sum(all_diff_restr([3 5 6 7],dumm_area & ismember(all_mk,'M'))<=thr_sig)~=0)]';
    sat_units_X(ar,:) = [sum(all_diff_restr([3],dumm_area & ismember(all_mk,'X'))<=thr_sig  & sum(all_diff_restr([5 6 7],dumm_area & ismember(all_mk,'X'))<=thr_sig)==0 ) ; ...
        sum(all_diff_restr([5],dumm_area & ismember(all_mk,'X'))<=thr_sig  & sum(all_diff_restr([6 7],dumm_area & ismember(all_mk,'X'))<=thr_sig)==0 ) ; ...
        sum(all_diff_restr([6],dumm_area & ismember(all_mk,'X'))<=thr_sig  & sum(all_diff_restr([5 7],dumm_area & ismember(all_mk,'X'))<=thr_sig)==0 ) ; ...
        sum(all_diff_restr([7],dumm_area & ismember(all_mk,'X'))<=thr_sig) ; ...
        sum(sum(all_diff_restr([3 5 6 7],dumm_area & ismember(all_mk,'X'))<=thr_sig)~=0)]';
end

sum(sat_units)
for ar = 1 : length(area2test)
[tbl,chi2_mot(ar),p_mot(ar)] = chi2_fms(sat_units(ar,1),sat_units(ar,end),sat_units(ar,end)-sat_units(ar,1),sat_units(ar,end))
end

test = {'SAT' 'SATxPB' 'SATxJU' 'SATxPBxJU' };

xax_ticks = ((0:1/(length(area2test)-1):1)-.5)*.55;
figure;
b = bar((sat_units(:,1:end-1)./repmat(sat_units(:,end),1,4))');hold on
for ar = 1 : length(b)
    b(ar).FaceColor = colorsArea(ar,:)
end
% for t = 1 : length(test)
%     plot(t+xax_ticks,rand_motif(:,t),'Color',[.6 .6 .6])
% end 
for t = 1 : length(test)
    bar_M = sat_units_M(:,t)./sat_units_M(:,end);
    bar_X = sat_units_X(:,t)./sat_units_X(:,end);
    plot(t+xax_ticks,bar_M,'o','MarkerSize',10,'MarkerFaceColor',[.8 .8 .8],'MarkerEdgeColor','k')
    plot(t+xax_ticks,bar_X,'v','MarkerSize',10,'MarkerFaceColor',[.65 .65 .65],'MarkerEdgeColor','k')
end
ylabel('Percent of SAT neurons')
set(gca,'XTick',1:length(test),'XTickLabel',test,'FontSize',16)

%- stats for common motif:
clear sat_units sat_units_M sat_units_X
all_diff_restr = all_diff(:,keep);

modeldata = table((all_diff_restr([3],:)<=thr_sig  & sum(all_diff_restr([5 6 7],:)<=thr_sig)==0)',...
                  (all_diff_restr([5],:)<=thr_sig  & sum(all_diff_restr([6 7],:)<=thr_sig)==0)',...
(all_diff_restr([6],:)<=thr_sig  & sum(all_diff_restr([5 7],:)<=thr_sig)==0)',...
(all_diff_restr([7],:)<=thr_sig)',all_area,all_mk,all_sess_id,'VariableNames',{'pref' 'pref_pb' 'pref_id' 'pref_pb_id' 'area' 'mk' 'sess'})


models_form = {'pref ~ 1 + area  + (1|mk) + (1|sess)' ; 'pref ~ 1 + area  + (1|mk)'};
%[lme,model_final] = model_comparison(modeldata,models_form,true);
lme = fitglme(modeldata,models_form{1}); model_final = models_form{1};

[pval,wald,thr_corr,pval_adj] = area_posthoc(lme,area2test,'y');
disp(['%%%%%%%%%% % sig ' test{t} ' with model = ' model_final ' %%%%%%%%%%'])
disp(anova(lme));disp(pval_adj);disp(wald);
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

models_form = {'pref_pb ~ 1 + area  + (1|mk) + (1|sess)' ; 'pref_pb ~ 1 + area  + (1|mk)'};
%[lme,model_final] = model_comparison(modeldata,models_form,true);
lme = fitglme(modeldata,models_form{1}); model_final = models_form{1};

[pval,wald,thr_corr,pval_adj] = area_posthoc(lme,area2test,'y');
disp(['%%%%%%%%%% % sig ' test{t} ' with model = ' model_final ' %%%%%%%%%%'])
disp(anova(lme));disp(pval_adj);disp(wald);
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

models_form = {'pref_id ~ 1 + area  + (1|mk) + (1|sess)' ; 'pref_id ~ 1 + area  + (1|mk)'};
%[lme,model_final] = model_comparison(modeldata,models_form,true);
lme = fitglme(modeldata,models_form{1}); model_final = models_form{1};

[pval,wald,thr_corr,pval_adj] = area_posthoc(lme,area2test,'y');
disp(['%%%%%%%%%% % sig ' test{t} ' with model = ' model_final ' %%%%%%%%%%'])
disp(anova(lme));disp(pval_adj);disp(wald);
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

models_form = {'pref_pb_id ~ 1 + area  + (1|mk) + (1|sess)' ; 'pref_pb_id ~ 1 + area  + (1|mk)'};
%[lme,model_final] = model_comparison(modeldata,models_form,true);
lme = fitglme(modeldata,models_form{1}); model_final = models_form{1};

[pval,wald,thr_corr,pval_adj] = area_posthoc(lme,area2test,'y');
disp(['%%%%%%%%%% % sig ' test{t} ' with model = ' model_final ' %%%%%%%%%%'])
disp(anova(lme));disp(pval_adj);disp(wald);
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')


%% plot examples PSTH - Figure 7A

%- for saving some examples..
% pathFig = 'C:\Users\Fred\Dropbox\Rudebeck Lab\Posters and Talks\PSTH Pref Example Neurons\';

%- 2 examples
s = [73 76];
n = [9 22];

for x = 1 : length(s)
    clearvars -except s n res_anova pathFig x
    if sum(res_anova(s(x)).Ap([3 5 6 7],n(x))<0.001)~=0

        %- extract the name of the file where the timestamps are!
        if res_anova(s(x)).neurons_info(n(x),1)<10; addzeros = '00';
        elseif res_anova(s(x)).neurons_info(n(x),1)<100; addzeros = '0';
        else addzeros = '';
        end
        if res_anova(s(x)).neurons_info(n(x),2)<10; addzeros2 = '0';
        else addzeros2 = '';
        end

        file_id = [res_anova(s(x)).session '_Ch' addzeros num2str(res_anova(s(x)).neurons_info(n(x),1)) ...
            '_Clus' addzeros2 num2str(res_anova(s(x)).neurons_info(n(x),2)) '.mat'];

        %- find the file in the different folder
        if strcmp(res_anova(s(x)).session(1),'M')
            path4file = utils_POTT_SPKfolder('MORBIER');
        elseif strcmp(res_anova(s(x)).session(1),'X') & datenum(res_anova(s(x)).session(2:7),'mmddyy')<datenum('090920','mmddyy')
            path4file = utils_POTT_SPKfolder('MIMIC1');
        else
            path4file = utils_POTT_SPKfolder('MIMIC2');
        end

        pref_ch = res_anova(s(x)).pref_ch;

        %- load unit
        Unit = load([path4file file_id]);
        Unit.info.area{1} = res_anova(s(x)).neurons_area{n(x)}; %- update area label based on histo

        %- load the event file from that session
        evt = load([path4file res_anova(s(x)).session '_behav.mat']);

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

        fighandle1 = rasterPSTH_2bins(Unit,evt2take,evtName,cond',2000,pre,post);
        %  fighandle1 =   rasterPSTH_2bins_SUB(Unit,evt2take,evtName,cond',2000,pre,post);
        pause
    end
end


