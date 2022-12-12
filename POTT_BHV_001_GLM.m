%% Main script for BEHAV analysis of the POTT dataset

%- to do:
%- update behav files with t_FPon... Can't compute the initiation time otherwise  // ADDED, just have to re-run all
%- update behav files mimic with Dynamic task juice info (for prevJuice) // DONE
%- check if dispersion flag should be true for main model, it is now for sliding window one.....
% clear
% path2go = utils_POTT_SPKfolder('MIMIC2')
% path2copy = 'C:\Users\Fred\Dropbox\Rudebeck Lab\ANA-POTT-BehavPrefChange\data\'
% 
% list = dir([path2go 'X*a_behav.mat']);
% for i = 1 : length(list)
%     disp(i)
%     copyfile([path2go list(i).name],path2copy)
% end


clear

mk = 'Mimic';

%- locate the files
path2go = ('C:\Users\Fred\Dropbox\Rudebeck Lab\ANA-POTT-BehavPrefChange\data\')
%path2go = ('/Users/Fred/Dropbox/Rudebeck Lab/ANA-POTT-BehavPrefChange/data/')

%- list all the files
if strcmp(mk(1:2),'Mi')
    list = dir([path2go 'X*a_behav.mat']);
else
    list = dir([path2go 'M*a_behav.mat']);
end

perm = false;
nPerm = 50;
showfig = false;

%- reorganize list by date
days = [];
for ss = 1 : length(list)
    days = [days;datenum(list(ss).name(2:7),'mmddyy') , ss];
end
date_order = sortrows(days);
list = list(date_order(:,2));

for sess = 1 : length(list) % for every sessions
    filename = [path2go list(sess).name];

    [ALL(sess,1),param] = behav_model(filename,false,showfig); %- no permutations
    
    if perm
        for p = 1 : nPerm
            [ALLperm(sess,p),~] = behav_model(filename,perm,showfig);
        end
    end
    
end
if perm
    % save([path2go mk '_behav_norm_DEVAL_prevJuice_ft.mat'],'ALL','ALLperm','param','-v7.3')
    save([path2go mk '_behav_norm_ALL_prevJuice_ft.mat'],'ALL','ALLperm','param','-v7.3')
else
    save([path2go mk '_behav_norm_ALL_prevJuice_ft.mat'],'ALL','param','-v7.3')
end


%% POTT - Behavioral analyses and Figure 1

clear
cd('C:\Users\Fred\Dropbox\Rudebeck Lab\ANA-POTT-BehavPrefChange\data\')
%cd('/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-BehavPrefChange/data/')
%cd('\\10.81.115.14\fred\POTT\data\behav')

M1 = load('Morbier_behav_norm_ALL_prevJuice_only.mat','ALL','param')
M2 = load('Mimic_behav_norm_ALL_prevJuice_only.mat','ALL','param')

thr_trend = 0.05;

%- remove sessions before threshold for Morbier, if needed
subset = {'052318' '010122'} ; %- for anything after last change
if ~isempty(subset)
    lim = datenum(subset,'mmddyy');
    for ss = 1 : length(M1.ALL)
        days(ss) = datenum(M1.ALL(ss).name(end-16:end-11),'mmddyy');
    end
    takeme = (days>=lim(1) & days<lim(2)) ;
else
    takeme = true(1,length(M1.ALL));
end
M1.ALL = M1.ALL(takeme);


%- save a matrix with the names of the considered sessions + converge model
clear behav_sess
x=0;
for m = 1 : 2
    clear conv
    eval(['ALL = M' num2str(m) '.ALL;']);
   
    takeme = [ALL(:).converge]==1;
    sessions = cat(1,ALL(takeme).name);
    for s = 1 : size(sessions,1)
        x = x + 1;
        behav_sess{x,1} = sessions(s,end-17:end-10);
    end
end

save('Sessions.mat','behav_sess')


%% Fig 1B - Behav model R squared + Estimates
figure;
colors = cbrewer('qual', 'Paired', 10);

subplot(1,6,1)
%- plot adjusted R2
for m = 1 : 2
    clear converge r_adj  keep
    eval(['ALL = M' num2str(m) '.ALL;']);
    for i = 1 : length(ALL)
        r_adj(:,i) = ALL(i).mdl.Rsquared.Adjusted;
        converge(:,i) = ALL(i).converge;
    end
    keep = converge==1 ;
    if m == 1
        col = colors(3:4,:);
        mm = -0.2;
    else
        col = colors(9:10,:);
        mm = .2;
    end
    X= r_adj(1,keep);
    yl=1+mm;
    wdth = .35;
    boxplot_ind(X,yl,wdth,col)
end

set(gca,'view',[90 -90],'color','none','FontSize',16);
set(gca,'YTick',[])
xlabel('R-adj')

xlim([0 1])

subplot(1,6,3:6)
%- plot the estimates
for m = 1 : 2
    clear converge tStat pVal Estimate r_adj LogLik nTr keep
    eval(['ALL = M' num2str(m) '.ALL;']);
    for i = 1 : length(ALL)
        tStat(:,i) = ALL(i).mdl.Coefficients.tStat;
        pVal(:,i) = ALL(i).mdl.Coefficients.pValue;
        Estimate(:,i) = ALL(i).mdl.Coefficients.Estimate;
        r_adj(:,i) = ALL(i).mdl.Rsquared.Adjusted;
        LogLik(:,i) = ALL(i).mdl.LogLikelihood;
        converge(:,i) = ALL(i).converge;
        nTr(:,i) = height(ALL(i).T);
    end
    
    keep = converge==1  ;
    
    if m == 1
        col = colors(3:4,:);
        mm = -0.2;
        line([0 0],[0 6],'Color','k');hold on
    else
        col = colors(9:10,:);
        mm = .2;
    end
    for i  = 1 : length( ALL(1).mdl.CoefficientNames)-1
        X= Estimate(i+1,keep);
        yl=i+mm;
        wdth = .35;
        boxplot_ind(X,yl,wdth,col)
    end
    set(gca,'view',[90 -90],'color','none','FontSize',16);
    set(gca,'YTick',1:length(ALL(1).mdl.CoefficientNames)-1,'YTickLabel',ALL(1).mdl.CoefficientNames(2:end),'YTickLabelRotation',25)
    disp([sum(keep) length(keep)])
    if m == 1
        text(19.5,4,['mk M'],'FontSize',16,'Color',col(2,:))
    else
        text(18,4,['mk X'],'FontSize',16,'Color',col(2,:))
    end
end
xlabel('Estimates')
xlim([-20 20])
ylim([0 length( ALL(1).mdl.CoefficientNames)])



%% look at choice consistency in same juice trials for the different preferences
both_pref = [];
figure;
for m = 1 : 2
    subplot(1,2,m)
    clear converge tStat pVal Estimate r_adj LogLik nTr keep
    eval(['ALL = M' num2str(m) '.ALL;']);
    pref = [ALL(:).pref];
    p_trend = [ALL(:).p_trend];
    p_trend_bins = [ALL(:).p_trend_bins];
    choice_consistency_J1 = [ALL(:).choice_consistency_J1];
    choice_consistency_J2 = [ALL(:).choice_consistency_J2];
    converge = [ALL(:).converge];

    plot(pref,choice_consistency_J1,'o'); hold on
    plot(pref,choice_consistency_J2,'or');
    histogram(pref,[0:0.15:1],'Normalization','Probability')
    xlim([0 1])

    %- proportion of sessions with trend in residuals
    disp(['Monkey ' num2str(m)])
    disp('     nSig      nTot     %sig')
    disp([sum(p_trend(converge==1)<thr_trend) length(p_trend(converge==1)) 100*(sum(p_trend(converge==1)<thr_trend)/length(p_trend(converge==1)))])

    both_pref = [both_pref , pref];
end



