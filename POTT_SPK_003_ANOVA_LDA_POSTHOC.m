%% POST HOC FOR ANOVA + LDA = FIGURE 2 and Suppl. Figs

%- Author: Fred M. Stoll, Icahn School of Medicine at Mount Sinai, NY
%- Date: 2023.03

clear

path2go = '/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-BehavPrefChange/data/neurons/subset-final/'

%- load the anova results
% load([path2go 'res_ANOVA_150tr_full.mat'])
load([path2go 'res_ANOVA_full.mat'])
% load([path2go 'res_ANOVA_150tr_v3.mat']) %- without chosen side
params = [3 1];

%- load the decoding results
measures = {'I_chosenjuice' 'I_chosenproba'}% 'chosenproba'}
LDA = load([path2go 'res_LDA_stim.mat'])


thr=0.01;nsig=3; %- threshold for sig

%- for plotting over time
step = 20; % step in time for plotting purposes only (write down value every X time bin)
timesel = [2 4 5 6];
timesel_sub = [-.4 .98 ; 0 .98 ; -.2 .38 ; -.2 .8];
areas = utils_POTT_areas;
%area2test = {'12r' '12m' '12o' '12l' 'a11ml' '13l' '13m' 'AMG' 'LAI' '44' '45' 'Striatum'};
 area2test = {'vlPFC' 'OFC' 'IFG' 'LAI' 'AMG' };

%- COLOR ASSIGNMENT
order = [3 1 2 5 4 7];
colorsArea = cbrewer('qual', 'Set2', 8);
colorsArea = colorsArea(order,:);
colorsArea_sub = cbrewer('qual', 'Pastel2', 8);
colorsArea_sub = colorsArea_sub(order,:);

if length(area2test)>6
    colorsArea =     [210 221 244 ; 162 184 232 ; 107 142 218 ; 41 91 202 ;
        %191 238 223 ; 125 221 181 ; 69 207 163 ; 8 191 133 ;
        185 247 193 ; 70 179 144 ; 0 142 98 ;
        % 185 247 193 ; 141 217 193 ; 70 179 144 ; 0 142 98 ;
        231 138 195 ;
        166 216 84 ;
        255 187 114 ; 255 102 43 ;
        229 196 148]/255
end

%- periods for computing average number of sig neurons during time bins : name / time alignment / times
% periods = {'REF'  , 'Stim_onset' , [-600 -100]  ;
periods = {'REF'  , 'FixFP_onset' , [100 700]  ;
    'STIM' , 'Stim_onset' ,  [200 800] ;
    'REW'  ,  'Rew_onset' ,  [100 700] };
%- second column refers to times_evts matrix, like: times_evts = {'FixFP_onset' 'Stim_onset' 'Resp_onset' 'FixResp_onset' 'FB_onset' 'Rew_onset' 'FB_offset'};

all_diff = [];
all_diff_omega = [];
all_diff_betas = [];
all_sess = [];
all_units = [];
all_mk = [];
all_neuron_id = [];
for s = 1 : length(res_anova)
    if ~isempty(res_anova(s).neurons_area)
        all_diff_omega = [all_diff_omega ; res_anova(s).ins_diff.Omega2];
        all_diff_betas = [all_diff_betas ; res_anova(s).ins_diff.Betas];
        all_diff = [all_diff ; res_anova(s).ins_diff.Ap];
        all_units = [all_units ; res_anova(s).neurons_area];
        all_mk = [all_mk ; repmat({res_anova(s).session(1)},size(res_anova(s).neurons_area,1),1)];
        all_sess = [all_sess ; repmat(s,size(res_anova(s).neurons_area,1),1)];

        for n = 1 : length(res_anova(s).neurons_area)
            %- extract the name of the file where the timestamps are!
            if res_anova(s).neurons_info(n,1)<10; addzeros = '00';
            elseif res_anova(s).neurons_info(n,1)<100; addzeros = '0';
            else addzeros = '';
            end
            if res_anova(s).neurons_info(n,2)<10; addzeros2 = '0';
            else addzeros2 = '';
            end

            all_neuron_id = [all_neuron_id ; [res_anova(s).session '_Ch' addzeros num2str(res_anova(s).neurons_info(n,1)) ...
                '_Clus' addzeros2 num2str(res_anova(s).neurons_info(n,2)) ]];
        end

    end
end

for p = 1 : size(all_diff,2)
    dumm = squeeze(all_diff(:,p,:));
    pval_sig = zeros(size(dumm));
    for k = 1 : length(dumm(:,1)) %- for each neuron, check if below thr
        [~,idxs] = findenough(dumm(k,:),thr,nsig,'<=');
        pval_sig(k,idxs)=1;clear idxs
    end
    all_diff_sig(:,p,:) = pval_sig;
end


%- remove neurons where models fails
keep = (sum(isnan(squeeze(all_diff(:,1,:)))')==0)';
mk1 = ismember(all_mk,'M')

all_units_area = cell(size(all_units));
for ar = 1 : length(area2test)
    eval(['takeit = ismember(all_units,areas.' area2test{ar} ');'])
    all_units_area(takeit) = area2test(ar);
end

%- main figure (Fig 2)

params = [3 1];


%-
%- number of units // sessions (WHEN ENOUGH TRIALS FOR FIG2!!!)
 for ar = 1 : length(area2test)
     eval(['takeit = ismember(all_units,areas.' area2test{ar} ') & keep;'])
     nbunit_tot(ar,:) = [sum(takeit & mk1) sum(takeit & ~mk1)];
  %   nbsess_tot(ar,:) = [length(unique(all_sess(takeit & mk1))) length(unique(all_sess(takeit & ~mk1)))];
 end

figure;m = 0;

for p = 1 : length(params)
    m = m+1;
    
    % barplot percent sig
    for bins = 1:2 % : size(periods,1)
        time_chunk = bins_considered == find(ismember(times_evts,periods{bins,2}));
        time_considered = find(time_chunk & time>=periods{bins,3}(1) & time<=periods{bins,3}(2));
        for ar = 1 : length(area2test)
            eval(['takeit = ismember(all_units,areas.' area2test{ar} ') & keep;'])
            dumm = squeeze(all_diff_sig(:,params(p),:))';
            
            bar_both(ar,:) = [sum(sum(dumm(time_considered,takeit))~=0) sum(takeit)];
            bar_M(ar,:) = [sum(sum(dumm(time_considered,takeit & mk1 ))~=0) sum(takeit & mk1)];
            bar_X(ar,:) = [sum(sum(dumm(time_considered,takeit & ~mk1 ))~=0) sum(takeit & ~mk1)];

            if bins == 2
                name_sig{p} = all_neuron_id(sum(dumm(time_considered,:))~=0 & keep',:);
            end
        end
        if bins == 1
            ref_sig = bar_both(:,1)./bar_both(:,2);
        end
        
        %- extract percent sig and stat
        if bins == 2
            sig = (sum(squeeze(all_diff_sig(keep,params(p),time_considered))')~=0)';
            
            modeldata_su = table(sig,all_units_area(keep),all_mk(keep),categorical(all_sess(keep)), 'VariableNames',{'sig' 'area' 'mk' 'sess'});
            models_form = {'sig ~ 1 + area  + (1|mk) + (1|sess)' ; 'sig ~ 1 + area  + (1|mk)'};
            % [lme,model_final] = model_comparison(modeldata_su,models_form,true);
            lme = fitglme(modeldata_su,models_form{1},'Distribution','Binomial'); model_final = models_form{1};

            [~,wald_percNeurons,~,pval_adj_percNeurons]  = area_posthoc(lme,area2test,'n');
            disp(['%%%%%%%%%% Percent of ' res_anova(1).ins_diff.predic{params(p)} ' neurons with model = ' model_final ' %%%%%%%%%%'])
            disp(anova(lme));disp(pval_adj_percNeurons);disp(wald_percNeurons);
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        end
        
    end
    
    subplot(2,5,m)
    for ar = 1 : length(area2test)
        bar(ar, [bar_both(ar,1)./bar_both(ar,2)]*100,'FaceColor',colorsArea(ar,:));hold on
    end
    plot(ref_sig*100,'Color',[.6 .6 .6])
    plot([bar_M(:,1)./bar_M(:,2)]*100,'o','MarkerSize',10,'MarkerFaceColor',[.8 .8 .8],'MarkerEdgeColor','k')
    plot([bar_X(:,1)./bar_X(:,2)]*100,'v','MarkerSize',10,'MarkerFaceColor',[.65 .65 .65],'MarkerEdgeColor','k')
    set(gca,'Xtick',1:length(area2test),'XtickLabel',area2test,'XtickLabelRotation',30,'FontSize',16)
    ylim([0 80])
    xlim([0 length(area2test)+1])
    ylabel('Percent significant neurons')
    
    %- time resolved omega2
    m = m + 1;
    subplot(2,5,[m m+1])
    modeldata_o2=table(); %- for stats
    for ar = 1 : length(area2test)
        eval(['takeit = ismember(all_units,areas.' area2test{ar} ') & keep;'])
        
        %- show perc sig
        % dumm = squeeze(all_diff_sig(:,params(p),:))';
        % perc = mean(dumm(:,takeit),2)*100;
        %- show Exp Variance in sig neurons
        time_chunk = bins_considered == find(ismember(times_evts,periods{2,2}));
        time_considered = find(time_chunk & time>=periods{2,3}(1) & time<=periods{2,3}(2));
        sig_anytime = ismember(bins_considered,timesel);
        
        sig = squeeze(all_diff_sig(:,params(p),:))';
        dumm = squeeze(all_diff_omega(:,params(p),:))';
        perc = nanmean(dumm(:,sum(sig(time_considered,:))>0 & takeit'),2);
        perc_sem = nanstd(dumm(:,sum(sig(time_considered,:))>0 & takeit')')' / sqrt(sum(takeit));

        om2(p,ar,1) = mean(mean(dumm(time_considered,sum(sig(time_considered,:))>0 & takeit'))   );
        om2(p,ar,2) = std(mean(dumm(time_considered,sum(sig(time_considered,:))>0 & takeit'))   )/sqrt(sum(takeit) );
        om2_m(p,ar) = mean(mean(dumm(time_considered,sum(sig(time_considered,:))>0 & takeit' & mk1'  ))   );
        om2_x(p,ar) = mean(mean(dumm(time_considered,sum(sig(time_considered,:))>0 & takeit' & ~mk1'  ))   );
       
        cd_curr = sum(sig(time_considered,:))>0 & takeit';
        modeldata_o2 = [modeldata_o2 ; table(double(mean(dumm(time_considered,cd_curr)))',all_units_area(cd_curr),all_mk(cd_curr),categorical(all_sess(cd_curr)), 'VariableNames',{'omega2' 'area' 'mk' 'sess'})];

        %- extract the Beta values
        sig1 = squeeze(all_diff_sig(:,params(1),:))';
        sig2 = squeeze(all_diff_sig(:,params(2),:))';
        dumm = squeeze(all_diff_betas(:,params(p),:))';
        dumm(~sig1 & ~sig2)=NaN;
        betas{p,ar} = nanmean(dumm(time_considered,takeit));
        betas_sigN{ar} = [sum(sig1(time_considered,takeit))>0 ; sum(sig2(time_considered,takeit))>0 ;...
                        sum(sig1(time_considered,takeit))>0 & sum(sig2(time_considered,takeit))>0 ];

%         V = dumm(time_considered,takeit);
%         [A,X] = max(abs(V));
%         betas{p,ar} = [];
%         for i = 1 : length(A)
%             betas{p,ar} = [betas{p,ar} A(i).*sign(V(X(i),i))];
%         end

        V = dumm(time_considered,sum(sig1(time_considered,:))>0 & sum(sig2(time_considered,:))>0 & takeit');
        betas_sig{p,ar} = mean(V);
%         [A,X] = max(abs(V));
%         betas_sig{p,ar} = [];
%         for i = 1 : length(A)
%             betas_sig{p,ar} = [betas_sig{p,ar} A(i).*sign(V(X(i),i))];
%         end
        
        gaps = [1 find(diff(time)~=subsp)+1 length(time)];
        for t = 1 : length(timesel)
            time_considered = find(bins_considered == timesel(t) & time>=1000*timesel_sub(t,1) & time<=1000*timesel_sub(t,2));
            plot_me{t} = time_considered;
        end
        timestart = 1;
        xlab = [];
        for t = 1 : length(timesel)
            timeend = timestart+length(plot_me{t})-1;
            plot(timestart:timeend,perc(plot_me{t}),'Color',colorsArea(ar,:),'LineWidth',2); hold on
            plot(timestart:timeend,perc(plot_me{t})+perc_sem(plot_me{t}),'Color',colorsArea(ar,:),'LineWidth',.5); hold on
            plot(timestart:timeend,perc(plot_me{t})-perc_sem(plot_me{t}),'Color',colorsArea(ar,:),'LineWidth',.5); hold on
            line([timeend timeend],[0 80],'Color',[.6 .6 .6])
            timestart = timeend;
            xlab = [xlab time(plot_me{t})];
        end
        per = find(time(plot_me{1})>=periods{2,3}(1) & time(plot_me{1})<=periods{2,3}(2));
        line([per(1) per(end)],[0 0],'Color','k','LineWidth',4)
           
        set(gca,'Xtick',1:step:length(xlab),'XtickLabel',xlab(1:step:end)/1000,'XtickLabelRotation',30,'FontSize',16)
        text(10,.1-(ar/200),[area2test{ar} ' - ' num2str(sum(takeit)) ],'Color',colorsArea(ar,:),'FontSize',16)
        
    end
    title(res_anova(1).ins_diff.predic{params(p)})
    ylabel('Mean Omega2')
    xlabel('Time (sec)')
    xlim([0 timeend+1])
    ylim([0 .1])
 
    %- stats
    models_form = {'omega2 ~ 1 + area  + (1|mk) + (1|sess)' ; 'omega2 ~ 1 + area  + (1|mk)'};
    %[lme,model_final] = model_comparison(modeldata_o2,models_form,false);
    lme = fitglme(modeldata_o2,models_form{1}); model_final = models_form{1};

    [~,wald_omegaNeurons,~,pval_adj_omegaNeurons]  = area_posthoc(lme,area2test,'n');
    disp(['%%%%%%%%%% Variance of ' res_anova(1).ins_diff.predic{params(p)} ' neurons with model = ' model_final ' %%%%%%%%%%'])
    disp(anova(lme));disp(pval_adj_omegaNeurons);disp(wald_omegaNeurons);
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')


    if length(area2test)<=6
        
        %- decoding
        m = m +2;
        subplot(2,5,m);hold on
        
        dumm_areas = [];
        maxPerf = [];
        sessions_all=char();
        for ar = 1 : length(area2test)
            perf_all=[];
            nUnit=[];
            sessions=char();
            
            eval(['curr = LDA.res_LDA.' measures{p} '.' area2test{ar} ';']);
            for s = 1 : length(curr)
                nUnit(s) = length(curr(s).takeit);
                perf_all(s,:)=nanmean(curr(s).perf,2);
                sessions(s,:) = curr(s).lda_sess(1:7);
            end
            
            %- decoding only done on bin 2 = so time to take [200-800]
            %t_lim = length(time_considered);
            
            nUnit_all{p,ar} = nUnit;
            %  maxPerf = [maxPerf ; [nanmean(perf_all(:,end-t_lim:end)')' , ar*ones(length(nUnit),1) , nUnit']];
            maxPerf = [maxPerf ; [perf_all , ar*ones(length(nUnit),1) , nUnit']];
            sessions_all = [sessions_all ; sessions];
            mk = sessions_all(:,1);
            mk = ismember(mk,'M')+1;
            dumm_areas = [dumm_areas ; repmat(area2test(ar),length(nUnit),1)];
        end
                
        %- plot decoding perf
        mm = 0;
        for i  = 1 : length( area2test)
            X = maxPerf(maxPerf(:,2)==i,1);
            mk_sub = mk(maxPerf(:,2)==i);
            yl=i+mm;
            wdth = .5;
         %   boxplot_ind_mk(X,yl,wdth,[.8 .8 .8 ; .65 .65 .65 ; colorsArea(i,:)],mk_sub)
               boxplot_ind(X,yl,wdth,[colorsArea_sub(i,:) ; colorsArea(i,:)])
               nbsess(i,:) = [sum(mk_sub==2) sum(mk_sub==1)];

        end
        set(gca,'view',[90 -90],'color','none','FontSize',16);
        set(gca,'YTick',1:length(area2test),'YTickLabel',area2test,'YTickLabelRotation',30)
        ylim([0 length(area2test)+1]);
        title(measures{p})
        xlabel('Decoding probability')
        
        %- significance
        modeldata_pop = table(maxPerf(:,1),dumm_areas,maxPerf(:,3),sessions_all,sessions_all(:,1),'VariableNames',{'perf' 'area' 'nb' 'sess' 'mk'});
        models_form = {'perf ~ 1 + area  + (1|mk) + (1|sess)' ; 'perf ~ 1 + area  + (1|mk)'};
        % [lme,model_final] = model_comparison(modeldata_pop,models_form,false);
        lme = fitglme(modeldata_pop,models_form{1}); model_final = models_form{1};

        [pval,wald,thr_corr,pval_adj] = area_posthoc(lme,area2test,'n');
        disp(['%%%%%%%%%% Decoding perf with model = ' model_final ' %%%%%%%%%%'])
        disp(anova(lme));disp(pval_adj);disp(wald);
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        
        %- plot the significance
        pval2=[pval_adj , zeros(length(area2test),1)]+ [pval_adj' ; zeros(1,length(area2test))];
        for ar1  = 1 : length( area2test)
            Xmax = max(maxPerf(maxPerf(:,2)==ar1,1));
            Xmin = min(maxPerf(maxPerf(:,2)==ar1,1));
            Xmean = mean(maxPerf(maxPerf(:,2)==ar1,1));
            updown=[1.5 2];
            for ar2 = 1 : length(area2test)
                Xmean2 = mean(maxPerf(maxPerf(:,2)==ar2,1));
                if ar1~=ar2 & pval2(ar1,ar2)<thr_corr & Xmean<Xmean2
                    text(Xmax+(updown(1)/120),ar1,'*','Color',colorsArea(ar2,:),'FontSize',20,'FontWeight','bold','HorizontalAlignment','center')
                    updown(1) = updown(1) + 1;
                elseif ar1~=ar2 & pval2(ar1,ar2)<thr_corr & Xmean>Xmean2
                    text(Xmin-(updown(2)/120),ar1,'*','Color',colorsArea(ar2,:),'FontSize',20,'FontWeight','bold','HorizontalAlignment','center')
                    updown(2) = updown(2) + 1;
                end
            end
        end
        
        %- sig decoding proportion
        m = m + 1;
        subplot(2,5,m);hold on
        
        dumm_perf=[];dumm_mk =[];dumm_sess=[];dumm_area=[];
        for ar = 1 : length(area2test)
            eval(['curr = LDA.res_LDA.' measures{p} '.' area2test{ar} ';']);
            perf = [curr(:).perf_pval];
            bar_decod(ar,:) = [sum(perf<0.05) length(perf<0.05) ];
            
            sess = cat(1,curr(:).lda_sess);
            mk = sess(:,1);
            
            bar_decodM(ar,:) = [sum(perf(mk=='M')<0.05) length(perf(mk=='M')) ];
            bar_decodX(ar,:) = [sum(perf(mk=='X')<0.05) length(perf(mk=='X')) ];
            dumm_perf = [dumm_perf ; perf'<0.05 ];
            dumm_mk = [dumm_mk ; mk ];
            dumm_sess = [dumm_sess ; sess ];
            dumm_area = [dumm_area ; repmat(area2test(ar),size(mk)) ];
        end
        for ar = 1 : length(area2test)
            bar(ar, [bar_decod(ar,1)./bar_decod(ar,2)]*100,'FaceColor',colorsArea(ar,:));hold on
        end
        plot([bar_decodM(:,1)./bar_decodM(:,2)]*100,'o','MarkerSize',10,'MarkerFaceColor',[.8 .8 .8],'MarkerEdgeColor','k')
        plot([bar_decodX(:,1)./bar_decodX(:,2)]*100,'v','MarkerSize',10,'MarkerFaceColor',[.65 .65 .65],'MarkerEdgeColor','k')
        set(gca,'Xtick',1:length(area2test),'XtickLabel',area2test,'XtickLabelRotation',30,'FontSize',16)
        ylim([0 100])
        xlim([0 length(area2test)+1])
        ylabel('Percent significant decoding')
        
        %- significance percent sig decoding
            
        modeldata_pop = table(dumm_perf,dumm_area,dumm_mk,dumm_sess, 'VariableNames',{'sig' 'area' 'mk' 'sess'});
        models_form = {'sig ~ 1 + area  + (1|mk) + (1|sess)' ; 'sig ~ 1 + area  + (1|mk)'};
        % [lme,model_final] = model_comparison(modeldata_pop,models_form,true);
        lme = fitglme(modeldata_pop,models_form{1},'Distribution','Binomial'); model_final = models_form{1};

        [pval,wald_popsig_adj,thr_corr,pval_popsig_adj] = area_posthoc(lme,area2test,'n');
        disp(['%%%%%%%%%% Percent sig decoding with model = ' model_final ' %%%%%%%%%%'])
        disp(anova(lme));disp(pval_popsig_adj);disp(wald_popsig_adj);
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        
        %- Build Table s2
        tab_s2_part1 =[];ar_s2 =[];
        for ar1 = 1 : length(area2test)-1
            for ar2 = 1 : length(area2test)
                if ar2>ar1
                    n=1;pp=0;
                    while pp == 0
                        n = n + 1;
                        pp = round(pval_adj_percNeurons(ar2,ar1),n);
                    end
                    pp = round(pval_adj_percNeurons(ar2,ar1),n+1);
                   
                    tab_s2_part1 = [tab_s2_part1;{num2str(round(wald_percNeurons(ar2,ar1),2)) num2str(pp)}];
                    ar_s2 = [ar_s2 ;{area2test{ar1} area2test{ar2}}];
                end
            end
        end
        tab_s2_part2 =[];ar_s2 =[];
        for ar1 = 1 : length(area2test)-1
            for ar2 = 1 : length(area2test)
                if ar2>ar1
                    n=1;pp=0;
                    while pp == 0
                        n = n + 1;
                        pp = round(pval_adj(ar2,ar1),n);
                    end
                    pp = round(pval_adj(ar2,ar1),n+1);
                   
                    tab_s2_part2 = [tab_s2_part2;{num2str(round(wald(ar2,ar1),2)) num2str(pp)}];
                    ar_s2 = [ar_s2 ;{area2test{ar1} area2test{ar2}}];
                end
            end
        end
        tab_s2_part3 =[];ar_s3 =[];
        for ar1 = 1 : length(area2test)-1
            for ar2 = 1 : length(area2test)
                if ar2>ar1
                    n=1;pp=0;
                    while pp == 0
                        n = n + 1;
                        pp = round(pval_popsig_adj(ar2,ar1),n);
                    end
                    pp = round(pval_popsig_adj(ar2,ar1),n+1);
                   
                    tab_s2_part3 = [tab_s2_part3;{num2str(round(wald_popsig_adj(ar2,ar1),2)) num2str(pp)}];
                    ar_s3 = [ar_s3 ;{area2test{ar1} area2test{ar2}}];
                end
            end
        end
        if p == 1
            tab_s2 = [tab_s2_part1 tab_s2_part2 tab_s2_part3];
        else
            tab_s2 = [tab_s2 tab_s2_part1 tab_s2_part2 tab_s2_part3];
        end
        
    else
        m = m + 3;

    end
    

end
set(gcf, 'Color', [1 1 1]);

% save('C:\Users\Fred\Dropbox\Rudebeck Lab\ANA-POTT-BehavPrefChange\POTT_sigUnits_name.mat','name_sig')

grpstats(modeldata_pop,{'area' 'mk'},"numel")

figure;
for m = 1 : 2
    subplot(1,2,m)
    for ar = 1 : length(area2test)
        bar(ar, squeeze(om2(m,ar,1)) ,'FaceColor',colorsArea(ar,:));hold on
        line([ar ar], [squeeze(om2(m,ar,1))-squeeze(om2(m,ar,2)) squeeze(om2(m,ar,1))+squeeze(om2(m,ar,2))] ,'Color',colorsArea(ar,:));hold on
        plot(ar,om2_m(m,ar),'o','MarkerSize',10,'MarkerFaceColor',[.8 .8 .8],'MarkerEdgeColor','k')
        plot(ar,om2_x(m,ar),'v','MarkerSize',10,'MarkerFaceColor',[.65 .65 .65],'MarkerEdgeColor','k')
    end
    ylim([0 0.07])    
    set(gca,'Xtick',1:length(area2test),'XtickLabel',area2test,'XtickLabelRotation',30,'FontSize',16)
    xlim([0 length(area2test)+1])
    ylabel('Mean Omega2')

end



% figure;
% for ar = 1 : length(area2test)
%     subplot(1,length(area2test),ar)
%     plot(betas{2,ar}, betas{1,ar} ,'.', 'Color',colorsArea(ar,:));hold on
%     plot(betas_sig{2,ar}, betas_sig{1,ar} ,'o', 'Color',colorsArea(ar,:));hold on
%     xlim([-.25 .25]);
%     ylim([-.25 .25])
%     [r,p]=corrcoef(betas{2,ar}, betas{1,ar})
% end

figure;
for ar = 1 : length(area2test)
    subplot(3,length(area2test),ar)
    plot(betas{2,ar}(betas_sigN{ar}(1,:)==1), betas{1,ar}(betas_sigN{ar}(1,:)==1) ,'.', 'Color',colorsArea(ar,:));hold on
    xlim([-.25 .25]);ylim([-.25 .25])
    subplot(3,length(area2test),ar+length(area2test))
    plot(betas{2,ar}(betas_sigN{ar}(2,:)==1), betas{1,ar}(betas_sigN{ar}(2,:)==1) ,'.', 'Color',colorsArea(ar,:));hold on
    xlim([-.25 .25]);ylim([-.25 .25])
    subplot(3,length(area2test),ar+length(area2test)+length(area2test))
    plot(betas{2,ar}(betas_sigN{ar}(3,:)==1), betas{1,ar}(betas_sigN{ar}(3,:)==1) ,'.', 'Color',colorsArea(ar,:));hold on
    xlim([-.25 .25]);ylim([-.25 .25])
end
%     plot(betas_sig{2,ar}, betas_sig{1,ar} ,'o', 'Color',colorsArea(ar,:));hold on
%     [r,p]=corrcoef(betas{2,ar}, betas{1,ar})
% end

  sig1 = squeeze(all_diff_sig(:,params(1),:))';
  sig2 = squeeze(all_diff_sig(:,params(2),:))';
  figure;
   for ar = 1 : length(area2test)
       subplot(1,length(area2test),ar)
            eval(['takeit = ismember(all_units,areas.' area2test{ar} ') & keep;'])
      ju = mean(sig1(:,takeit)==1 & sig2(:,takeit)==0,2) ;
      both = mean(sig1(:,takeit)==1 & sig2(:,takeit)==1,2) ;
      pb = mean(sig1(:,takeit)==0 & sig2(:,takeit)==1,2) ;
plot(pb(bins_considered==2));hold on;
plot(ju(bins_considered==2));
plot(both(bins_considered==2));
plot(both(bins_considered==2)+ju(bins_considered==2)+pb(bins_considered==2));
ylim([0 .5])
   end



%        sig1 = squeeze(all_diff_sig(:,params(1),:))';
%         sig2 = squeeze(all_diff_sig(:,params(2),:))'; 
%         dumm1 = squeeze(all_diff_betas(:,params(1),:))';
%         dumm2 = squeeze(all_diff_betas(:,params(2),:))';
%         om1 = squeeze(all_diff_omega(:,params(1),:))';
%         om2 = squeeze(all_diff_omega(:,params(2),:))';
%         dumm(~sig1 & ~sig2)=NaN;
% 
% 
% figure;plot(abs(om2(150,:)),abs(dumm2(150,:)),'.')

%%
figure;
time_chunk = bins_considered == find(ismember(times_evts,periods{2,2}));
time_considered = find(time_chunk & time>=0 & time<=1000);
subplot(2,1,1)
line([0 40],[0 40],'Color','k');hold on
for ar = 1 : length(area2test)
    eval(['takeit = ismember(all_units,areas.' area2test{ar} ') & keep ;'])
 %   eval(['takeit = ismember(all_units,areas.' area2test{ar} ') & keep & ismember(all_mk,"M");'])
    
    dumm = squeeze(all_diff_sig(:,params(1),:))';
    dumm2 = squeeze(all_diff_sig(:,params(2),:))';
    for t = 1 : length(time_considered)
        perc_ju(t) = sum(dumm(time_considered(t),takeit))/sum(takeit);
        perc_pb(t) = sum(dumm2(time_considered(t),takeit))/sum(takeit);
    end
    
    last = max([find(perc_ju==max(perc_ju),1,'last') find(perc_pb==max(perc_pb),1,'last')]);
    
    plot(perc_pb(1:last)*100,perc_ju(1:last)*100,'.-','Color',colorsArea(ar,:),'LineWidth',2,'MarkerSize',10)
    hold on
    tt = time(time_considered)==300
    plot(perc_pb(tt)*100,perc_ju(tt)*100,'.','Color',colorsArea(ar,:),'MarkerSize',30)
end
box on
xlabel('Percent Proba neurons')
ylabel('Percent Juice neurons')
grid on
set(gca,'FontSize',16);
xlim([0 40])
ylim([0 40])
    
subplot(2,1,2)
line([0 .08],[0 .08],'Color','k');hold on
for ar = 1 : length(area2test)
    eval(['takeit = ismember(all_units,areas.' area2test{ar} ') & keep ;'])
%    eval(['takeit = ismember(all_units,areas.' area2test{ar} ') & keep & ismember(all_mk,"M");'])

    sig = squeeze(all_diff_sig(:,params(1),:))';
    sig2 = squeeze(all_diff_sig(:,params(2),:))';
    dumm = squeeze(all_diff_omega(:,params(1),:))';
    dumm2 = squeeze(all_diff_omega(:,params(2),:))';
    for t = 1 : length(time_considered)
        perc_ju(t) = nanmean(dumm(time_considered(t),sum(sig(time_considered,:))>0 & sum(sig2(time_considered,:))>0 & takeit'));
        perc_pb(t) = nanmean(dumm2(time_considered(t),sum(sig(time_considered,:))>0 & sum(sig2(time_considered,:))>0 & takeit'));
    end
    
    last = max([find(perc_ju==max(perc_ju),1,'last') find(perc_pb==max(perc_pb),1,'last')]);
    
    plot(perc_pb(1:last),perc_ju(1:last),'.-','Color',colorsArea(ar,:),'LineWidth',2,'MarkerSize',10)
    hold on
    tt = time(time_considered)==300
    plot(perc_pb(tt),perc_ju(tt),'.','Color',colorsArea(ar,:),'MarkerSize',30)
end
box on
xlabel('Mean Omega Proba')
ylabel('Mean Omega Juice')
grid on
set(gca,'FontSize',16);
xlim([0 .08])
ylim([0 .08])
          title('Proba + Juice neurons')



%% 

time_chunk = bins_considered == find(ismember(times_evts,periods{2,2}));
time_considered = find(time_chunk & time>=0 & time<=1000);

sig_time=[];
sig_time_peak=[];
for p = 1 : length(params)
    dumm = squeeze(all_diff_sig(:,params(p),:))';
    dumm_omega = squeeze(all_diff_omega(:,params(p),:))';
    for n = 1 : size(dumm,2)
        [tt,ttt] = findenough(dumm(time_considered,n)',1,nsig,'==');
        if isempty(tt)
            sig_time(n,p)=NaN;
            sig_time_peak(n,p)=NaN;
        else
            sig_time(n,p)=tt(1);
            tt2 = find(dumm_omega(time_considered(ttt),n)==max(dumm_omega(time_considered(ttt),n)),1,'first');
            sig_time_peak(n,p)=ttt(tt2(1));
        end
    end
end


%% add single monkey points...


time_chunk = bins_considered == find(ismember(times_evts,periods{2,2}));
time_considered = find(time_chunk & time>=periods{2,3}(1) & time<=periods{2,3}(2));
time_considered = find(time_chunk & time>=0 & time<=1000);
%- some ref during FP onset
% time_chunk = session.bins_considered == find(ismember(session.times_evts,{'FixFP_onset'}));
% time_considered = find(time_chunk & session.time>=periods{2,3}(1) & session.time<=periods{2,3}(2));

    clear ratio_perm

for ar = 1 : length(area2test)
    eval(['take_me = ismember(all_units,areas.' area2test{ar} ') & keep & ismember(all_mk,"M");'])
  %  eval(['take_me = ismember(all_units,areas.' area2test{ar} ') & keep ;'])
    dumm_juice = squeeze(all_diff_sig(:,params(1),:))';
    dumm_proba = squeeze(all_diff_sig(:,params(2),:))';
    sig_proba = sum(dumm_proba(time_considered,take_me))~=0;
    sig_juice = sum(dumm_juice(time_considered,take_me))~=0;
    
    vennlike = [sum(sig_juice(1,:)==1 & sig_proba(1,:)==0)  ...
                sum(sig_juice(1,:)==1 & sig_proba(1,:)==1) ...
                sum(sig_juice(1,:)==0 & sig_proba(1,:)==1) ...
                sum(sig_juice(1,:)==0 & sig_proba(1,:)==0)];
    
    prop_sig(ar,:) = vennlike/sum(vennlike);
%     OR_tab = [sum(sig_juice(1,:)==0 & sig_proba(1,:)==0) sum(sig_juice(1,:)==0 & sig_proba(1,:)==1) ;
%               sum(sig_juice(1,:)==1 & sig_proba(1,:)==0) sum(sig_juice(1,:)==1 & sig_proba(1,:)==1) ];
%     OR(ar) =    ( OR_tab(1,1) * OR_tab(2,2) ) / ( OR_tab(2,1) * OR_tab(1,2) );
    
    ratio(ar) = mean(sig_proba(sig_juice(1,:)==1)) / mean(sig_proba(sig_juice(1,:)==0)); %- P(Pb|Ju) / P(Pb|¬Ju)
  %  ratio(ar) = mean(sig_juice(sig_proba(1,:)==1)) / mean(sig_juice(sig_proba(1,:)==0)); %- P(Pb|Ju) / P(Pb|¬Ju)
    ratio(ar) = (ratio(ar)-1)*100 ; %- put that in percent change

    %- chi2 method
   [tbl,chi2stat(ar),pvals(ar)] = chi2_fms(vennlike(2),sum(vennlike(1:2)),vennlike(3),sum(vennlike(3:4)));

%     ratio2(ar) = mean(sig_juice(sig_proba(1,:)==1)) / mean(sig_juice(sig_proba(1,:)==0)); %- P(Ju|Pb) / P(Ju|¬Pb)
%     ratio2(ar) = (ratio2(ar)-1)*100 ; %- put that in percent change
%     
%    [tbl,chi2stat2(ar),pvals2(ar)] = chi2_fms(vennlike(2),sum(vennlike(2:3)),vennlike(1),sum(vennlike([1 4])));

   %- permutation method
    nP = 1000;
    for p = 1 : nP
        permuted = randperm(length(sig_proba));
        permAll = [sig_juice(1,:)' sig_proba(1,permuted)'];

        ratio_perm(ar,p) = mean(permAll(permAll(:,1)==1,2)) / mean(permAll(permAll(:,1)==0,2)) ; %- P(Pb|Ju) / P(Pb|¬Ju)
        ratio_perm(ar,p) = (ratio_perm(ar,p)-1)*100 ; %- put that in percent change

    end
    
    
    
end

figure;
subplot(3,10,[7:10]);hold on
for ar = 1 : length(area2test)
     b = barh(ar,ratio(ar));
     text(ratio(ar)+3,ar,['p = ' num2str(pvals(ar))],'FontSize',14)
     b.FaceColor = colorsArea(ar,:);
     plot(ratio_perm(ar,:),ar+(0.7*rand(length(ratio_perm(ar,:)),1))-.35,'.','Color',[.6 .6 .6])
end
xlim([-50 200])  ;box on;axis ij
set(gca,'Ytick',1:length(area2test),'YtickLabel',area2test,'FontSize',16)    
xlabel('Percent Change = P(Pb|Ju) / P(Pb|¬Ju)')    
    
subplot(3,10,[1:5]);hold on
colorsBar = [cbrewer('qual', 'Accent', 3) ; .6 .6 .6];
b = barh(prop_sig,'stacked');
box on;axis ij
set(gca,'Ytick',1:length(area2test),'YtickLabel',area2test,'FontSize',16)    
xlabel('Proportion of neurons')    
xpos = [.05 .25 .55 .75]
xlab = {'Juice' 'Juice+Proba' 'Proba' 'NS'}
for i = 1 : 4
    b(i).FaceColor = colorsBar(i,:);
    text(xpos(i),.25,xlab{i},'Color',colorsBar(i,:),'FontSize',16)
end

%- latency of proba or identity
sig_tab = [];

    lat_dumm = sig_time(keep,:);
    sess_dumm = all_sess(keep,:);
    %perc(perc<time_considered(1) | perc >time_considered(end))=NaN;
    mk = mk1(keep);
    ar_dumm = all_units_area(keep);
    %- take the one that are sig for both Juice and Proba
 %   perc_cosig = lat_dumm(~isnan(lat_dumm(:,1)) & ~isnan(lat_dumm(:,2)),[1 2]);
 %   mk_cosig = mk(~isnan(lat_dumm(:,1)) & ~isnan(lat_dumm(:,2))  )  ;
 %   sess_cosig = sess_dumm(~isnan(lat_dumm(:,1)) & ~isnan(lat_dumm(:,2))  )  ;
    subtime = time(time_considered);

    lat_dumm_ms = NaN(size(lat_dumm));
    for i = 1 : length(lat_dumm)
        if ~isnan(lat_dumm(i,1))
            lat_dumm_ms(i,1) = [subtime(lat_dumm(i,1))];
        end
        if ~isnan(lat_dumm(i,2))
            lat_dumm_ms(i,2) = [subtime(lat_dumm(i,2))];
        end        
    end
    modeldata_all = table(lat_dumm_ms(:,1),lat_dumm_ms(:,2),categorical(mk),categorical(sess_dumm),ar_dumm,'VariableNames',{'lat_ju','lat_pb','mk','sess','area'});

        models_form = {'lat_pb ~ 1 + area  + (1|mk) + (1|sess)' ; 'lat_pb ~ 1 + area  + (1|mk)'};
        % [lme,model_final] = model_comparison(cosig_tab,models_form,false);
        lme_pb = fitglme(modeldata_all(~isnan(modeldata_all.lat_pb),:),models_form{1}); model_final = models_form{1};

        [pval,wald_latpb_adj,thr_corr,pval_latpb_adj] = area_posthoc(lme_pb,area2test,'y');
        disp(['%%%%%%%%%% Latency sig pb with model = ' model_final ' %%%%%%%%%%'])
        disp(anova(lme_pb));disp(pval_latpb_adj);disp(wald_latpb_adj);
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        
        emm = emmeans(lme_pb,{'area'});
        emm.table

        models_form = {'lat_ju ~ 1 + area  + (1|mk) + (1|sess)' ; 'lat_ju ~ 1 + area  + (1|mk)'};
        % [lme,model_final] = model_comparison(cosig_tab,models_form,false);
        lme_ju = fitglme(modeldata_all(~isnan(modeldata_all.lat_ju),:),models_form{1}); model_final = models_form{1};

        [pval,wald_latju_adj,thr_corr,pval_latju_adj] = area_posthoc(lme_ju,area2test,'y');
        disp(['%%%%%%%%%% Latency sig ju with model = ' model_final ' %%%%%%%%%%'])
        disp(anova(lme_ju));disp(pval_latju_adj);disp(wald_latju_adj);
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        
        emm = emmeans(lme_ju,{'area'});
        emm.table


%% NOT USED
%- cosig difference in latencies
cosig_tab = [];

%figure;
for ar = 1 : length(area2test)
    eval(['take_me = ismember(all_units,areas.' area2test{ar} ') & keep ;'])
    perc = sig_time(take_me,:);
    sess_dumm = all_sess(take_me,:);
    %perc(perc<time_considered(1) | perc >time_considered(end))=NaN;
    mk = mk1(take_me);
   
    %- take the one that are sig for both Juice and Proba
    perc_cosig = perc(~isnan(perc(:,1)) & ~isnan(perc(:,2)),[1 2]);
    mk_cosig = mk(~isnan(perc(:,1)) & ~isnan(perc(:,2))  )  ;
    sess_cosig = sess_dumm(~isnan(perc(:,1)) & ~isnan(perc(:,2))  )  ;
    subtime = time(time_considered);
    x = perc_cosig(:,2);
    y = perc_cosig(:,1);
    for i = 1 : length(x)
        x(i,:) = [subtime(x(i,1))];
        y(i,:) = [subtime(y(i,1))];
    end
    
    cosig_tab = [cosig_tab ; table(x-y,mk_cosig,sess_cosig,repmat(area2test(ar),size(mk_cosig)),'VariableNames',{'lat','mk','sess','area'})];

    %     [r,p] = corrcoef(x,y)
    % allRP(ar,:) = [r(2,1) , p(2,1)];
    
    xgrid=linspace(0,600,100);
    ygrid=xgrid;
    [x1,y1] = meshgrid(xgrid, ygrid);
    % Perform kernel density estimate
    % [x y] is actual data, xi is the desired grid points to evaluate
    % f is an estimate of the density, ep(:,1) is the X location, ep(:,2) is the y location
    xi = [x1(:) y1(:)];
    [f,ep]=ksdensity([x y],xi); % remove the outputs to see a 3D plot of the distribution
    % format data in matrix for contourf and plot
    X = reshape(ep(:,1),length(xgrid),length(ygrid));
    Y = reshape(ep(:,2),length(xgrid),length(ygrid));
    Z = reshape(f,length(xgrid),length(ygrid));
    %  subplot(1,length(area2test),ar);
    %contourf(X,Y,Z,10)
    %  imagesc(X(1,:),Y(:,1),Z);axis xy
    
    % N = size(Z(:),1);                                      % Number of ‘Experiments’ In Data Set
    % yMean = mean(Z(:));                                    % Mean Of All Experiments At Each Value Of ‘x’
    % ySEM = std(Z(:));                              % Compute ‘Standard Error Of The Mean’ Of All Experiments At Each Value Of ‘x’
    % CI95 = tinv([.95], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    % yCI95 = bsxfun(@times, ySEM, CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’
    % contour(X,Y,Z,[yCI95 yCI95],'Color',colorsArea(ar,:))
    [a,b] = hist(x-y,[-1000:20:1000]);
% % subplot(1,2,2);
% % 
% % plot(b,cumsum(a)/sum(a),'Color',colorsArea(ar,:),'LineWidth',2);
% % hold on
% % text(-500,1-(ar/30),[area2test{ar} ],'Color',colorsArea(ar,:),'FontSize',16)

subplot(3,10,[11:15 21:25]);
    
    contour(X,Y,Z,[0.6*max(Z(:)) 0.6*max(Z(:))],'Color',colorsArea(ar,:),'LineWidth',1)
    %contour(X,Y,Z,[0.7*max(Z(:)) 0.7*max(Z(:))],'Color',colorsArea(ar,:),'LineWidth',1)
    hold on
    %contour(X,Y,Z,[.995*max(Z(:)) .995*max(Z(:))],'Color',colorsArea(ar,:),'LineWidth',3)
    [xx,yy]=find(Z==max(Z(:)))
    plot(X(xx,yy),Y(xx,yy),'.','Color',colorsArea(ar,:),'MarkerSize',20)
text(50,600-(ar*20),[area2test{ar} ],'Color',colorsArea(ar,:),'FontSize',16)
%     hax = axes('Position', [.35, .35, .3, .3]);
% bar(hax,y,'EdgeColor','none')
% set(hax,'XTick',[])

  end
  
subplot(3,10,[11:15 21:25]);
line([0 600],[0 600],'Color','k')
xlabel('Latency of Chosen Probability encoding (ms)')
ylabel('Latency of Chosen Juice encoding (ms)')
set(gca,'FontSize',16)
title('Contour (60% of the max) of the bivariate distribution + max')

% % subplot(1,2,2);
% % line([-1000 1000],[.5 .5],'Color','k')
% % line([0 0],[0 1],'Color','k')
% % xlabel('Latency difference (Chosen Proba - Chosen Juice) (ms)')
% % ylabel('Cumulative proportion of neurons')
% % set(gca,'FontSize',16)
            
        models_form = {'lat ~ 1 + area  + (1|mk) + (1|sess)' ; 'lat ~ 1 + area  + (1|mk)'};
        % [lme,model_final] = model_comparison(cosig_tab,models_form,false);
        lme = fitglme(cosig_tab,models_form{1}); model_final = models_form{1};

        [pval,wald_lat_adj,thr_corr,pval_lat_adj] = area_posthoc(lme,area2test,'y');
        disp(['%%%%%%%%%% Latency co-sig with model = ' model_final ' %%%%%%%%%%'])
        disp(anova(lme));disp(pval_lat_adj);disp(wald_lat_adj);
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        
     
%lme_lat = fitglme(cosig_tab,'lat ~ 1 + area + (1|mk)')
%anova(lme_lat)
%[pval,wald_lat,thr,pval_adj_lat] = area_posthoc(lme_lat,area2test,'n')

emm = emmeans(lme,{'area'});
emm.table

%- some comparison to 0 (not used...)
clear pval_diff0 zval_diff0
for ar = 1 : length(area2test)
    [P,H,STATS] = signtest(cosig_tab.lat(ismember(cosig_tab.area,area2test{ar}) ));
    pval_diff0(ar,:) = P ;
    zval_diff0(ar,:) = STATS.zval;
end
[h, crit_p, adj_p]=fdr_bh(pval_diff0(:),.05);

%%

% cosig_xtab = [];
% cosig_ytab = [];
% 
% %figure;
% for ar = 1 : length(area2test)
%     eval(['take_me = ismember(all_units,areas.' area2test{ar} ') & keep;'])
%     perc = sig_time(take_me,:);
%     %perc(perc<time_considered(1) | perc >time_considered(end))=NaN;
%     mk = mk1(take_me);
%    
%     %- take the one that are sig for both Juice and Proba
%     perc_cosig = perc(~isnan(perc(:,1)) & ~isnan(perc(:,2)),[1 2]);
%     mk_cosig = mk(~isnan(perc(:,1)) & ~isnan(perc(:,2))  )  ;
%     subtime = time(time_considered);
%     x = perc_cosig(:,2);
%     y = perc_cosig(:,1);
%     for i = 1 : length(x)
%         x(i,:) = [subtime(x(i,1))];
%         y(i,:) = [subtime(y(i,1))];
%     end
%     
%     cosig_xtab = [cosig_xtab ; table(x,mk_cosig,repmat(area2test(ar),size(mk_cosig)),'VariableNames',{'lat','mk','area'})];
%     cosig_ytab = [cosig_ytab ; table(y,mk_cosig,repmat(area2test(ar),size(mk_cosig)),'VariableNames',{'lat','mk','area'})];
% end
% 
% lme_x = fitglme(cosig_xtab,'lat ~ 1 + area + (1|mk)')
% lme_y = fitglme(cosig_ytab,'lat ~ 1 + area + (1|mk)')
% anova(lme_x)
% anova(lme_y)
% [pval,wald_x,thr,pval_adj_x] = area_posthoc(lme_x,area2test,'y')
% [pval,wald_y,thr,pval_adj_y] = area_posthoc(lme_y,area2test,'y')

