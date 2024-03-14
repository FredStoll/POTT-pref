%% POTT Pref - Cross-subspace decoding post-processing
%-
%- Analyses for Fig 4C-D of Stoll & Rudebeck, Neuron, 2024
%-
%- Author: Fred M. Stoll, Icahn School of Medicine at Mount Sinai, NY
%- Date: 2023.12
%- Related to: Stoll & Rudebeck, Neuron, 2024

clear

path2go = '/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-BehavPrefChange/data/neurons/subset-final/'

%- load the decoding results
measures = {'I_juice_by_chosenproba'  };% 'chosenproba'}
name = {'perf_pb' 'perf_fl' }
load([path2go 'res_LDA_cross_revision2_final.mat']);

area2test = {'vlPFC' 'OFC' 'IFG' 'LAI' 'AMG'}
order = [3 1 2 5 4 7];
colorsArea = cbrewer('qual', 'Set2', 8);
colorsArea = colorsArea(order,:);
colorsArea_sub = cbrewer('qual', 'Pastel2', 8);
colorsArea_sub = colorsArea_sub(order,:);

util.statout = @(out,p2take_name) ['F(' num2str(out.DF1(strcmp(out.Term,p2take_name))) ',' num2str(out.DF2(strcmp(out.Term,p2take_name))) ...
    ')=' num2str(out.FStat(strcmp(out.Term,p2take_name))) ...
    ', p(' p2take_name ')=' num2str(out.pValue(strcmp(out.Term,p2take_name)))];

maxPerf = [];
sessions_all=char();
for m = 1 : length(measures)
    tab(m).modeldata=[];

    for ar = 1 : length(area2test)
        eval(['curr = res_LDA.' measures{m} '.' area2test{ar}])
        ss = size(curr(1).perf);
        perf = [];
        sess=[];
        mk = [];
        perf_pval=[];
        for s = 1 : length(curr)
            perf(s,:) = nanmean(curr(s).perf,length(ss))';
            perf_pval(s,:) = curr(s).perf_pval';

            sess{s,1} = curr(s).lda_sess ;
            mk{s,1} = curr(s).lda_sess(1) ;

        end
        perf_sig = [curr(:).perf_pval];
        perf_sig= perf_sig<0.05;
        curr_area = repmat(area2test(ar),size(mk));

        if size(perf,2)==1
            tab(m).modeldata = [tab(m).modeldata ; table(perf,perf_sig',perf_raw_sig',curr_area,sess,mk,'VariableNames',{'perf' 'sig' 'raw_sig' 'area' 'session' 'mk'})];
        else
            tab(m).modeldata = [tab(m).modeldata ; table(perf(:,1),perf(:,2),perf(:,3),perf(:,4),perf_pval(:,1),perf_pval(:,2),perf_pval(:,3),perf_pval(:,4),curr_area,sess,mk,'VariableNames',{'perf_ju_projju' 'perf_ju_projpb' 'perf_pb_projju' 'perf_pb_projpb' 'pval1' 'pval2' 'pval3' 'pval4' 'area' 'session' 'mk'})];
        end
    end
end

data = tab(m).modeldata;
figure
for ar = 1 : length(area2test)
    subplot(2,length(area2test),ar)
    line([0 1],[0 1],'Color','k')
    hold on
    plot(data.perf_ju_projju(ismember(data.area,area2test{ar})),data.perf_ju_projpb(ismember(data.area,area2test{ar})),'.','Color',colorsArea(ar,:),'MarkerSize',15)
    xlim([.3 .9]);ylim([.3 .9]);box on

    diff_perf(ar,1) = nanmedian(data.perf_ju_projju(ismember(data.area,area2test{ar}))-data.perf_ju_projpb(ismember(data.area,area2test{ar})));
    xlabel('Flavor Decoding - Proj Flavor');
    ylabel('Flavor Decoding - Proj Proba');
    title(['Avg Perf decrease of ' num2str(diff_perf(ar)*100) '%'])

    subplot(2,length(area2test),ar+length(area2test))
    line([0 1],[0 1],'Color','k')
    hold on
    plot(data.perf_pb_projju(ismember(data.area,area2test{ar})),data.perf_pb_projpb(ismember(data.area,area2test{ar})),'.','Color',colorsArea(ar,:),'MarkerSize',15)
    xlim([.2 .7]);ylim([.2 .7]);box on

    diff_perf(ar,1) = nanmedian(data.perf_pb_projju(ismember(data.area,area2test{ar}))-data.perf_pb_projpb(ismember(data.area,area2test{ar})));
    xlabel('Proba Decoding - Proj Flavor');
    ylabel('Proba Decoding - Proj Proba');
    title(['Avg Perf decrease of ' num2str(diff_perf(ar)*100) '%'])
end

data = tab(m).modeldata;
psig=[];
for ar = 1 : length(area2test)
    for ppp=1:4
        eval(['pv = data.pval' num2str(ppp) '(ismember(data.area,area2test{ar}));'])
        eval(['pv_M = data.pval' num2str(ppp) '(ismember(data.area,area2test{ar}) & ismember(data.mk,"M"));'])
        eval(['pv_X = data.pval' num2str(ppp) '(ismember(data.area,area2test{ar}) & ismember(data.mk,"X"));'])
        psig(ppp,ar) = sum(pv<0.05) / length(pv);
        psig_M(ppp,ar) = sum(pv_M<0.05) / length(pv_M);
        psig_X(ppp,ar) = sum(pv_X<0.05) / length(pv_X);
    end
end
figure;
subplot(1,2,1)
b = bar(psig(1:2,:)*100);hold on
for pp=1:2
    plot(pp + [-0.3:0.15:0.3],psig_M(pp,:)*100,'o','MarkerSize',10,'MarkerFaceColor',[.8 .8 .8],'MarkerEdgeColor','k')
    plot(pp + [-0.3:0.15:0.3],psig_X(pp,:)*100,'v','MarkerSize',10,'MarkerFaceColor',[.65 .65 .65],'MarkerEdgeColor','k')
end
for ar = 1 : length(area2test)
    b(ar).FaceColor = colorsArea(ar,:);
    set(gca,"XTick",1:4,"XTickLabel",{'Flavor subspace' 'Proba subspace' })
end
ylabel('Percent significant decoding')
legend([area2test 'mk M' 'mk X'])

% stat juice / psig difference
modeldata2 = data;
modeldata2.psig = modeldata2.pval1<0.05;
modeldata2.psig2 = modeldata2.pval2<0.05;

modeldata3 = table([modeldata2.psig; modeldata2.psig2],[modeldata2.area; modeldata2.area],[modeldata2.session; modeldata2.session],[modeldata2.mk ;modeldata2.mk],[ones(size(modeldata2.mk)) ;2*ones(size(modeldata2.mk))],'VariableNames',{'sig' 'area' 'session' 'mk' 'subspace'});
modeldata3.subspace=categorical(modeldata3.subspace);
models_form = {'sig ~ 1 + area + subspace + (1|mk) + (1|session)' };
lme = fitglme(modeldata3,models_form{1},'Distribution','Binomial'); model_final = models_form{1};

[pval,wald,thr_corr,pval_adj] = area_posthoc(lme,area2test,'n');
disp(anova(lme));disp(pval_adj);disp(wald);
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

title({'Chosen Flavor decoding',util.statout(anova(lme),'subspace')})

% same for proba decoding
subplot(1,2,2)
b = bar(psig(3:4,:)*100);hold on
for pp=1:2
    plot(pp + [-0.3:0.15:0.3],psig_M(pp+2,:)*100,'o','MarkerSize',10,'MarkerFaceColor',[.8 .8 .8],'MarkerEdgeColor','k')
    plot(pp + [-0.3:0.15:0.3],psig_X(pp+2,:)*100,'v','MarkerSize',10,'MarkerFaceColor',[.65 .65 .65],'MarkerEdgeColor','k')
end
for ar = 1 : length(area2test)
    b(ar).FaceColor = colorsArea(ar,:);
    set(gca,"XTick",1:2,"XTickLabel",{'Flavor subspace' 'Proba subspace' })
end
ylabel('Percent significant decoding')

%- stat for proba
modeldata2 = data;
modeldata2.psig = modeldata2.pval3<0.05;
modeldata2.psig2 = modeldata2.pval4<0.05;

modeldata3 = table([modeldata2.psig; modeldata2.psig2],[modeldata2.area; modeldata2.area],[modeldata2.session; modeldata2.session],[modeldata2.mk ;modeldata2.mk],[ones(size(modeldata2.mk)) ;2*ones(size(modeldata2.mk))],'VariableNames',{'sig' 'area' 'session' 'mk' 'subspace'});
modeldata3.subspace=categorical(modeldata3.subspace);
models_form = {'sig ~ 1 + area + subspace + (1|mk) + (1|session)' };
% [lme,model_final] = model_comparison(modeldata_pop,models_form,false);
lme = fitglme(modeldata3,models_form{1},'Distribution','Binomial'); model_final = models_form{1};

[pval,wald,thr_corr,pval_adj] = area_posthoc(lme,area2test,'n');
disp(anova(lme));disp(pval_adj);disp(wald);
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

title({'Chosen Proba decoding',util.statout(anova(lme),'subspace')})


%% more stats, but for area differences (+TABLE S4)

% stat juice / perf difference
modeldata2 = data;
modeldata2.perf = (modeldata2.perf_ju_projju-modeldata2.perf_ju_projpb)./modeldata2.perf_ju_projju;
models_form = {'perf ~ 1 + area  + (1|mk) + (1|session)' ; 'perf ~ 1 + area  + (1|mk)'};
% [lme,model_final] = model_comparison(modeldata_pop,models_form,false);
lme = fitglme(modeldata2,models_form{1}); model_final = models_form{1};

disp(['flavor: ' util.statout(anova(lme),'area')])

[pval,wald,thr_corr,pval_adj] = area_posthoc(lme,area2test,'y');
[~,table_wp1] = make_stattable(pval_adj,wald,area2test);

% stat proba / perf difference
modeldata2 = data;
modeldata2.perf = (modeldata2.perf_pb_projpb-modeldata2.perf_pb_projju)./modeldata2.perf_pb_projpb;
models_form = {'perf ~ 1 + area  + (1|mk) + (1|session)' ; 'perf ~ 1 + area  + (1|mk)'};
% [lme,model_final] = model_comparison(modeldata_pop,models_form,false);
lme = fitglme(modeldata2,models_form{1}); model_final = models_form{1};

disp(['proba: ' util.statout(anova(lme),'area')])

[pval,wald,thr_corr,pval_adj] = area_posthoc(lme,area2test,'y');
[ar_tab,table_wp2] = make_stattable(pval_adj,wald,area2test);

table_s4 = [ar_tab table_wp1 table_wp2];
disp(table_s4)

