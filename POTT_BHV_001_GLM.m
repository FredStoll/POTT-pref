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

mk = 'Mimic'; % 'Morbier'

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

%% show juice preference for both monkey (Fig 4A)
figure;
for m = 1 : 2
    clear converge tStat pVal Estimate r_adj LogLik nTr keep
    eval(['ALL = M' num2str(m) '.ALL;']);
    pref = [ALL(:).pref];
    if m == 1
        col = colors(3:4,:);
        mm = -0.2;
    else
        col = colors(9:10,:);
        mm = .2;
    end
    yl=1+mm;
    wdth = .35;
    boxplot_ind(pref,yl,wdth,col);hold on
    set(gca,'view',[90 -90],'color','none','FontSize',16);
    
end
xlim([0 1])
set(gca,'YTick',[])
xlabel('Juice Preference')


%% plot full correlation matrix (Fig 6C)
clear corrmat_both
for m = 1 : 2
    clear corrmat_all
    eval(['ALL = M' num2str(m) '.ALL;']);

    converge = [ALL(:).converge]==1;

    for d = 1 : size(ALL,1)
        corrmat_all(d,:,:)= ALL(d).corrmat;
    end
    dumm = squeeze(nanmean(abs(corrmat_all(converge,:,:))));

    if m == 1 ; corrmat_both(m,:,:) = dumm;
    elseif m ==2 ; corrmat_both(m,:,:) = dumm';
    end
end
corrmat_both = squeeze(nanmean(corrmat_both));

color_corr = cbrewer('seq', 'Greys', 64);
figure;
h = imagesc(corrmat_both,[0 0.6]);colormap(color_corr);colorbar;
for p1 = 1 : length(corrmat_both)
    for p2 = 1 : length(corrmat_both)
        if p1>p2
            text(p2,p1,['\it' num2str(round(100*corrmat_both(p1,p2))/100)],'HorizontalAlignment','center','Color',colors(4,:),'FontWeight','bold','FontSize',12)
        elseif p1<p2
            text(p2,p1,['\it' num2str(round(100*corrmat_both(p1,p2))/100)],'HorizontalAlignment','center','Color',colors(10,:),'FontWeight','bold','FontSize',12)
        end
    end
end
set(gca,'XTick',1:length(corrmat_both),'XTickLabel',M1.param.mat_name,'YTickLabelRotation',25)
set(gca,'YTick',1:length(corrmat_both),'YTickLabel',M1.param.mat_name,'XTickLabelRotation',25)
set(gca,'FontSize',14);
hold on;
line([length(corrmat_both)+1 0],[length(corrmat_both)+1 0],'Color',[0.6 0.6 0.6],'LineWidth',2)
h = pixelgrid(h);
h.Children(1).Color = [0.6 0.6 0.6];
%h.Children(2).Color = [31 120 180]/255;


%% Correlation btw Zval and preference - Fig 6B

%- https://en.wikipedia.org/wiki/Riemann_sum
%- not really the integral!
%- I divide by the length cos sessions widely different in length..
figure;
for m = 1:2
    subplot(2,1,m);
    line([0 0],[-15 15],'Color','k');
    hold on
    line([-0.5 0.5],[0 0],'Color','k'); box on;
    ref_length = 50;
    colors = cbrewer('qual', 'Paired', 10);

    clear devia_ref Z_trend p_trend converge
    eval(['ALL = M' num2str(m) '.ALL;']);

    converge = [ALL(:).converge]==1;
    p_trend = [ALL(:).p_trend];

    for d = 1 : size(ALL,1)
        Z_trend(d) = ALL(d).Z_trend;
        % devia_ref(d) = sum(ALL(d).pref_bins-mean(ALL(d).pref_bins(1:ref_length)))/length(ALL(d).pref_bins);
        % devia_ref(d) = trapz(ALL(d).pref_bins-mean(ALL(d).pref_bins(1:ref_length)),0:1/(length(ALL(d).pref_bins)-1):1);
        devia_ref(d) = mean(ALL(d).pref_bins-mean(ALL(d).pref_bins(1:ref_length)));
    end

    %- remove session that did not converge
    devia_ref=devia_ref(converge);
    p_trend=p_trend(converge);
    Z_trend=Z_trend(converge);

    if m == 1
        col= colors(3:4,:);
    else
        col= colors(9:10,:);
    end

    plot(devia_ref(p_trend>=thr_trend),Z_trend(p_trend>=thr_trend),'.','Color',col(1,:),'MarkerSize',15);hold on
    plot(devia_ref(p_trend<thr_trend),Z_trend(p_trend<thr_trend),'.','Color',col(2,:),'MarkerSize',15);hold on

    xx = devia_ref;
    yy = Z_trend;
    ww = ones(size(Z_trend)); %- not weighted!
    [R,p] = corrcoef(xx,yy);
    slope = R(2,1) * ( std(repelem(yy,ww)) / std(repelem(xx,ww)) );
    intercept = mean(repelem(yy,ww)) - ( slope * mean(repelem(xx,ww))); %- or that, the same thing : sum(W.*yy') - ( slope *  sum(W.*xx'));
    corrline = (xx * slope) + intercept;
    plot(xx,corrline,'-','LineWidth',2,'Color',col(2,:))

    %ylim([-12.5 12.5]);
    %xlim([-.35 .35]);
    ylabel('Z-value (Trend in residuals)')
    xlabel('Integral of Juice preference ')
    set(gca,'FontSize',12)
    text(-.325,-8.5,['R=' num2str(round(100*R(2,1))/100)],'FontSize',12,'Color',col(2,:))
    text(-.325,-10.5,['p=' num2str(p(2,1))],'FontSize',12,'Color',col(2,:))
end

%% Other correlation with RT and IT (Supp fig)

%- https://en.wikipedia.org/wiki/Riemann_sum
%- not really the integral! 
%- I divide by the length cos sessions widely different in length..
figure;
for m = 1:2
subplot(2,1,m);
line([0 0],[-15 15],'Color','k');
hold on
line([-0.5 0.5],[0 0],'Color','k'); box on;
ref_length = 50;
colors = cbrewer('qual', 'Paired', 10);

clear devia_ref Z_trend p_trend converge
    eval(['ALL = M' num2str(m) '.ALL;']);
    
    converge = [ALL(:).converge]==1;
    p_trend = [ALL(:).p_trend];
    
    for d = 1 : size(ALL,1)
        Z_trend(d) = ALL(d).Z_trend;
        % devia_ref(d) = sum(ALL(d).pref_bins-mean(ALL(d).pref_bins(1:ref_length)))/length(ALL(d).pref_bins);
        % devia_ref(d) = trapz(ALL(d).pref_bins-mean(ALL(d).pref_bins(1:ref_length)),0:1/(length(ALL(d).pref_bins)-1):1);
        devia_ref(d) = mean(ALL(d).ft_bins-mean(ALL(d).ft_bins(1:ref_length)));
    end

    %- remove session that did not converge
    devia_ref=devia_ref(converge);
    p_trend=p_trend(converge);
    Z_trend=Z_trend(converge);
    
    if m == 1
        col= colors(3:4,:);
    else
        col= colors(9:10,:);
    end
    
    plot(devia_ref(p_trend>=thr_trend),Z_trend(p_trend>=thr_trend),'.','Color',col(1,:),'MarkerSize',15);hold on
    plot(devia_ref(p_trend<thr_trend),Z_trend(p_trend<thr_trend),'.','Color',col(2,:),'MarkerSize',15);hold on
    
    xx = devia_ref;
    yy = Z_trend;
    ww = ones(size(Z_trend)); %- not weighted!
    [R,p] = corrcoef(xx,yy);
            slope = R(2,1) * ( std(repelem(yy,ww)) / std(repelem(xx,ww)) );
            intercept = mean(repelem(yy,ww)) - ( slope * mean(repelem(xx,ww))); %- or that, the same thing : sum(W.*yy') - ( slope *  sum(W.*xx'));
            corrline = (xx * slope) + intercept;
           plot(xx,corrline,'-','LineWidth',2,'Color',col(2,:))

%ylim([-12.5 12.5]);
xlim([-.2 .2]);
ylabel('Z-value (Trend in residuals)')
xlabel('Integral of IT ')
set(gca,'FontSize',12)
text(-.125,-8.5,['R=' num2str(round(100*R(2,1))/100)],'FontSize',12,'Color',col(2,:))
text(-.125,-10.5,['p=' num2str(p(2,1))],'FontSize',12,'Color',col(2,:))


end

%- https://en.wikipedia.org/wiki/Riemann_sum
%- not really the integral! 
%- I divide by the length cos sessions widely different in length..
figure;
for m = 1:2
subplot(2,1,m);
line([0 0],[-15 15],'Color','k');
hold on
line([-0.5 0.5],[0 0],'Color','k'); box on;
ref_length = 50;
colors = cbrewer('qual', 'Paired', 10);

clear devia_ref Z_trend p_trend converge
    eval(['ALL = M' num2str(m) '.ALL;']);
    
    converge = [ALL(:).converge]==1;
    p_trend = [ALL(:).p_trend];
    
    for d = 1 : size(ALL,1)
        Z_trend(d) = ALL(d).Z_trend;
        % devia_ref(d) = sum(ALL(d).pref_bins-mean(ALL(d).pref_bins(1:ref_length)))/length(ALL(d).pref_bins);
        % devia_ref(d) = trapz(ALL(d).pref_bins-mean(ALL(d).pref_bins(1:ref_length)),0:1/(length(ALL(d).pref_bins)-1):1);
        devia_ref(d) = mean(ALL(d).rt_bins-mean(ALL(d).rt_bins(1:ref_length)));
    end

    %- remove session that did not converge
    devia_ref=devia_ref(converge);
    p_trend=p_trend(converge);
    Z_trend=Z_trend(converge);
    
    if m == 1
        col= colors(3:4,:);
    else
        col= colors(9:10,:);
    end
    
    plot(devia_ref(p_trend>=thr_trend),Z_trend(p_trend>=thr_trend),'.','Color',col(1,:),'MarkerSize',15);hold on
    plot(devia_ref(p_trend<thr_trend),Z_trend(p_trend<thr_trend),'.','Color',col(2,:),'MarkerSize',15);hold on
    
    xx = devia_ref;
    yy = Z_trend;
    ww = ones(size(Z_trend)); %- not weighted!
    [R,p] = corrcoef(xx,yy);
            slope = R(2,1) * ( std(repelem(yy,ww)) / std(repelem(xx,ww)) );
            intercept = mean(repelem(yy,ww)) - ( slope * mean(repelem(xx,ww))); %- or that, the same thing : sum(W.*yy') - ( slope *  sum(W.*xx'));
            corrline = (xx * slope) + intercept;
           plot(xx,corrline,'-','LineWidth',2,'Color',col(2,:))

%ylim([-12.5 12.5]);
xlim([-.2 .2]);
ylabel('Z-value (Trend in residuals)')
xlabel('Integral of RT ')
set(gca,'FontSize',12)
text(-.125,-8.5,['R=' num2str(round(100*R(2,1))/100)],'FontSize',12,'Color',col(2,:))
text(-.125,-10.5,['p=' num2str(p(2,1))],'FontSize',12,'Color',col(2,:))


end

%- https://en.wikipedia.org/wiki/Riemann_sum
%- not really the integral! 
%- I divide by the length cos sessions widely different in length..
figure;
for m = 1:2
subplot(2,1,m);
line([0 0],[-15 15],'Color','k');
hold on
line([-0.5 0.5],[0 0],'Color','k'); box on;
ref_length = 50;
colors = cbrewer('qual', 'Paired', 10);

clear devia_ref Z_trend p_trend converge
    eval(['ALL = M' num2str(m) '.ALL;']);
    
    converge = [ALL(:).converge]==1;
    p_trend = [ALL(:).p_trend];
    
    for d = 1 : size(ALL,1)
        Z_trend(d) = mean(ALL(d).pref_bins-mean(ALL(d).pref_bins(1:ref_length)));
        % devia_ref(d) = sum(ALL(d).pref_bins-mean(ALL(d).pref_bins(1:ref_length)))/length(ALL(d).pref_bins);
        % devia_ref(d) = trapz(ALL(d).pref_bins-mean(ALL(d).pref_bins(1:ref_length)),0:1/(length(ALL(d).pref_bins)-1):1);
        devia_ref(d) = mean(ALL(d).ft_bins-mean(ALL(d).ft_bins(1:ref_length)));
    end

    %- remove session that did not converge
    devia_ref=devia_ref(converge);
    p_trend=p_trend(converge);
    Z_trend=Z_trend(converge);
    
    if m == 1
        col= colors(3:4,:);
    else
        col= colors(9:10,:);
    end
    
    plot(devia_ref(p_trend>=thr_trend),Z_trend(p_trend>=thr_trend),'.','Color',col(1,:),'MarkerSize',15);hold on
    plot(devia_ref(p_trend<thr_trend),Z_trend(p_trend<thr_trend),'.','Color',col(2,:),'MarkerSize',15);hold on
    
    xx = devia_ref;
    yy = Z_trend;
    ww = ones(size(Z_trend)); %- not weighted!
    [R,p] = corrcoef(xx,yy);
            slope = R(2,1) * ( std(repelem(yy,ww)) / std(repelem(xx,ww)) );
            intercept = mean(repelem(yy,ww)) - ( slope * mean(repelem(xx,ww))); %- or that, the same thing : sum(W.*yy') - ( slope *  sum(W.*xx'));
            corrline = (xx * slope) + intercept;
           plot(xx,corrline,'-','LineWidth',2,'Color',col(2,:))

%ylim([-12.5 12.5]);
ylim([-.35 .35]);
xlim([-.2 .2]);
ylabel('Integral of Juice Pref')
xlabel('Integral of IT ')
set(gca,'FontSize',12)
text(-.125,-0.2,['R=' num2str(round(100*R(2,1))/100)],'FontSize',12,'Color',col(2,:))
text(-.125,-0.3,['p=' num2str(p(2,1))],'FontSize',12,'Color',col(2,:))


end


%% Example session with change in preference (Mimic) - Fig. 6A

ALL = M2.ALL;

figure

d = 145; %145 15
% for d = 1 : 186
%     if ALL(d).p_trend<0.05
colors = cbrewer('div', 'RdBu', 128);
norm = @(data) -1+((data-min(data))*2)/(max(data)-min(data)) ;

probas = norm([10:1:90]);

clear sensitivity_bins pref_bins radj_bins ft_bins rt_bins mdl_bins resid
x = 0;
nb_bins = 3;
dumm = round(height(ALL(d).T)/nb_bins);
bins = [(1 : dumm : height(ALL(d).T)) height(ALL(d).T)]
for b = 1 : nb_bins
    disp(b)
    warning off
    mdl_bins = fitglm(ALL(d).T(bins(b):bins(b+1),:),M1.param.modelspec,'Distribution','binomial','Link','logit');
    ft = nanmedian(table2array(ALL(d).T(:,'ft')));
    rt = nanmedian(table2array(ALL(d).T(:,'rt')));

    for i = 1 : length(probas)
        tab = table(probas(i)*ones(length(probas),1),probas',nanmedian(ft)*ones(length(probas),1),nanmedian(rt)*ones(length(probas),1),-ones(length(probas),1),-ones(length(probas),1),'VariableNames',{'probJ1','probJ2','ft','rt','prevJuice','prevJuiceReceived'});
        tab.prevJuice = categorical(tab.prevJuice);
        tab.prevJuiceReceived = categorical(tab.prevJuiceReceived);
        if i == 1
            tab_all = tab;
        else
            tab_all = vertcat(tab_all,tab);
        end
    end
    newf_bins = [];
    try
        [newf_bins , ~] = predict(mdl_bins,tab_all);
        newf_bins = reshape(newf_bins,length(probas),length(probas));
    end
    x = x + 1;
    subplot(4,nb_bins,x)
    imagesc([10:1:90],[10:1:90],newf_bins);axis xy;colormap(colors)
end
subplot(4,nb_bins,[x+1 : 4*nb_bins])

bar(ALL(d).pref_bins,'FaceColor',[.8 .8 .8],'BaseValue',mean(ALL(d).pref_bins(1:ref_length)))
hold on
plot(ALL(d).pref_bins,'Color','k','LineWidth',1)

ylim([0 1]);xlim([0 length(ALL(d).pref_bins)])
title(num2str(d))
line([0 length(ALL(d).pref_bins)],[.5 .5],'Color','k')
line([0 length(ALL(d).pref_bins)],[mean(ALL(d).pref_bins(1:ref_length)) mean(ALL(d).pref_bins(1:ref_length))],'Color','r')
for i = 1 : length(bins)-2
    line([bins(i+1) bins(i+1)],[0 1],'Color',[.7 .7 .7])
end
xlabel('Trials')
ylabel('Juice Preference')
set(gca,'FontSize',14)
%plot((ALL(d).resid),'Color','k')
% pause
% hold off
%     end
% end

Z_trend(d)
devia_ref(d)

%% some extra stuff - not needed 


figure;
for m = 1 : 2
    clear corrmat_all p_trend
    eval(['ALL = M' num2str(m) '.ALL;']);
    for d = 1 : size(ALL,1)
        corrmat_all(d,:,:)= ALL(d).corrmat;
        p_trend(d) = ALL(d).p_trend;
    end
    res_corr = squeeze(corrmat_all(:,2,1));
    
    if m == 1
        col = colors(3:4,:);
        mm = -0.2;
        line([0 0],[0 6],'Color','k');hold on
    else
        col = colors(9:10,:);
        mm = .2;
    end
    
        ns = abs(res_corr(p_trend>thr_trend,1));
        sig = abs(res_corr(p_trend<thr_trend,1));
        [H, pValue(m), KSstat(m)] = kstest2(sig,ns);

        
        wdth = .35;
        boxplot_ind(ns(~isnan(ns)),1+mm,wdth,col)
        boxplot_ind(sig(~isnan(sig)),1+1+mm,wdth,col)
    set(gca,'view',[90 -90],'color','none','FontSize',16);
    set(gca,'YTick',1:2,'YTickLabel',{'ns' 'Trend<0.05'},'YTickLabelRotation',25)
    xlabel('Correlation Juice Pref and smoothed Residuals')
    if m == 1
        text(.1,1.25,['mk M - p(KS)=' num2str(round(100000*pValue(m))/100000)],'FontSize',16,'Color',col(2,:))
    else
        text(.05,1.25,['mk X - p(KS)=' num2str(round(100000*pValue(m))/100000)],'FontSize',16,'Color',col(2,:))
    end
end
ylim([0 3])

figure;
dd = [4 10];
for m = 1 : 2
    clear corrmat_all p_trend Z_trend
    eval(['ALL = M' num2str(m) '.ALL;']);
    for d = 1 : size(ALL,1)
        corrmat_all(d,:,:)= ALL(d).corrmat;
        p_trend(d) = ALL(d).p_trend;
        Z_trend(d) = ALL(d).Z_trend;
    end
    res_corr = squeeze(corrmat_all(:,2,1));
    [r,p]=corrcoef(abs(Z_trend),abs(res_corr))
    plot(abs(Z_trend),abs(res_corr),'.','Color',colors(dd(m),:)); hold on
end
xlabel('Z trend');ylabel('Corr Juice Pref vs Residuals')




