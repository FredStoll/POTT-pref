%% DECODING vs PREFERENCE (BETWEEN SESSION)

%- Author: Fred M. Stoll, Icahn School of Medicine at Mount Sinai, NY
%- Date: 2023.03

clear

clear
path2go = '/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-BehavPrefChange/data/neurons/subset-final/'; %- path where SPKpool files are!

measures = {'I_chosenjuice' 'I_chosenproba' 'I_unchosenproba'  'I_chosenside' 'I_chosenproba_juice' }
name = {'perf_ju' 'perf_pb' 'perf_unpb' 'perf_si' 'perf_pb1' 'perf_pb2' }
%measures = {'I_chosenside' 'I_unchosenproba' 'I_chosenproba_juice'}
%name = {'perf_ju' 'perf_pb' 'perf_pb1' 'perf_pb2'}
load([path2go 'res_LDA_stim.mat'])
area2test = {'vlPFC' 'OFC' 'IFG' 'LAI'}
order = [3 1 2 5 4 7];
colorsArea = cbrewer('qual', 'Set2', 8);
colorsArea = colorsArea(order,:);
colorsArea_sub = cbrewer('qual', 'Pastel2', 8);
colorsArea_sub = colorsArea_sub(order,:);
% figure;
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
        for s = 1 : length(curr)
            perf(s,:) = nanmean(curr(s).perf,length(ss))';

            sess{s,1} = curr(s).lda_sess ;
            mk{s,1} = curr(s).lda_sess(1) ;
        end
        curr_area = repmat(area2test(ar),size(mk));
        if size(perf,2)==1
            tab(m).modeldata = [tab(m).modeldata ; table(perf,curr_area,sess,mk,'VariableNames',{'perf' 'area' 'session' 'mk'})];
        else
            tab(m).modeldata = [tab(m).modeldata ; table(perf(:,1),perf(:,2),curr_area,sess,mk,'VariableNames',{'perf1' 'perf2' 'area' 'session' 'mk'})];
        end
    end
end

%- load behav
behav = load('/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-BehavPrefChange/data/Morbier_behav_bins.mat','ALL','param')
behav2 = load('/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-BehavPrefChange/data/Mimic_behav_bins.mat','ALL','param')
ALL = [behav.ALL ; behav2.ALL]
% pref_session = [ALL.pref];
clear days
pref_session = NaN(length(ALL),1);
trend_session = NaN(length(ALL),1);
for i = 1 : length(ALL)
    pref_session(i,1) = ALL(i).pref;
    trend_session(i,1) = ALL(i).Z_trend_bins;
    days{i,:} = [ALL(i).name(end-17:end-11) 'a'];
end


%- match decoding matrices and behav
all_sess = unique([tab(1).modeldata.session ; tab(2).modeldata.session ; tab(3).modeldata.session ; tab(4).modeldata.session ; tab(5).modeldata.session])  ;

x=0;

for i = 1 : length(all_sess)
    for ar = 1 : length(area2test)
        x=x+1;
        data(x,1:6)=NaN;
        pref_trend(x,1:2)=NaN;
        sess{x,1} = all_sess{i};
        mk{x,1} = all_sess{i}(1);
        curr_area{x,1} = area2test{ar};
        idx = [];
        for m = 1 : 5
            %  eval(['t' num2str(m) '= find(ismember(tab(m).modeldata.session,all_sess(i)) & ismember(tab(m).modeldata.area,area2test(ar)) );']);
            t = find(ismember(tab(m).modeldata.session,all_sess(i)) & ismember(tab(m).modeldata.area,area2test(ar)) );
            if ~isempty(t) & m~=5
                data(x,m)=tab(m).modeldata.perf(t);
            elseif ~isempty(t) & m==5
                data(x,m:m+1)=[tab(m).modeldata.perf1(t) tab(m).modeldata.perf2(t)];
            end
        end
        idx_behav = ismember(days,all_sess(i));
        pref_trend(x,:) = [pref_session(idx_behav) trend_session(idx_behav) ];
    end
end
modeldata = table(data(:,1),data(:,2),data(:,3),data(:,4),data(:,5),data(:,6),pref_trend(:,1),pref_trend(:,2),curr_area,sess,mk,'VariableNames',[name {'pref' 'trend' 'area' 'session' 'mk'}]);
rmv = sum(isnan(modeldata{:,1:6}),2)==6;
modeldata(rmv,:)=[];
data = modeldata;

%- reject AMG, not enough pref/neurons
modeldata(ismember(modeldata.area,'AMG'),:)=[];

%%
grp = {'modeldata.perf_ju' , 'abs(modeldata.pref-0.5)' ; ...
    'modeldata.perf_pb' , 'abs(modeldata.pref-0.5)' ;...
    'modeldata.perf_unpb' , 'abs(modeldata.pref-0.5)' ;...
    'modeldata.perf_si' , 'abs(modeldata.pref-0.5)' ;...
    'modeldata.perf_pb2-modeldata.perf_pb1' , 'modeldata.pref' }

colorsPr = cbrewer('seq', 'Greys', 100);
figure;
sub = 0;
for m = 1 : size(grp,1)
    y_all = eval([grp{m,1}]) ;
    x_all = eval([grp{m,2}]) ;

    ysc = [round(min(y_all),2)-.05 round(max(y_all),2)+.05];
    if strcmp(grp{m,2},'abs(modeldata.pref-0.5)')
        xsc = 0.5;
    else
        xsc = 1;
    end
    for ar = 1 : length(area2test)
        sub = sub + 1;
        subplot(size(grp,1),length(area2test),sub)
        take = ismember(modeldata.area,area2test(ar)) & ~isnan(x_all) & ~isnan(y_all);
        %   take = ismember(modeldata.area,area2test(ar)) & ismember(modeldata.mk,'X') & ~isnan(x_all) & ~isnan(y_all);
        x = x_all(take);
        y = y_all(take);

        [r,p] = corrcoef(x,y)
        allRP(ar,:) = [r(2,1) , p(2,1)];

        xgrid=linspace(0,xsc,100);
        ygrid=linspace(ysc(1),ysc(1),100);
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

        colormap(colorsPr)
        hold on
        line([0.5 0.5],[ysc(1) ysc(2)],'Color',[.6 .6 .6])
        line([0 1],[0 0],'Color',[.6 .6 .6])
        plot(x,y,'.','Color',colorsArea_sub(ar,:))
        text(.05,ysc(1)+0.02,['r(' num2str(length(x)-2) ')=' num2str(round(allRP(ar,1),3)) ' p=' num2str(round_pval(allRP(ar,2)))])

        cf = fit(x,y,'poly1'); % fit
        p_handle = plot(cf,'k','predfunc',.95); % Plot fit
        set(p_handle,'Color',colorsArea(ar,:),'LineWidth',2);
        legend off
        xlabel(grp{m,2});ylabel(grp{m,1})
        title(area2test{ar})
        xlim([0 xsc]); ylim([ysc(1) ysc(2)])
        box on


    end

end


abs(modeldata.pref-0.5)

tbl_juice = table(modeldata.perf_ju,abs(modeldata.pref-0.5),modeldata.area,modeldata.mk,'VariableNames',{'perf' 'pref' 'area' 'mk'})
lme_juice = fitglme(tbl_juice,'perf ~ 1 + area*pref + (1|mk)')

anova(lme_juice)
% lme_juice_fitted = fitted(lme_juice)
%
% figure;plot(lme_juice_fitted,modeldata.perf_ju,'o')


tbl_proba = table(modeldata.perf_pb,abs(modeldata.pref-0.5),modeldata.area,modeldata.mk,'VariableNames',{'perf' 'pref' 'area' 'mk'})
lme_proba = fitglme(tbl_proba,'perf ~ 1 + area*pref + (1|mk)')

anova(lme_proba)


tbl_pb12 = table(modeldata.perf_pb2-modeldata.perf_pb1,modeldata.pref,modeldata.area,modeldata.mk,'VariableNames',{'perf' 'pref' 'area' 'mk'})
lme_pb12 = fitglme(tbl_pb12,'perf ~ 1 + area*pref + (1|mk)')

anova(lme_pb12)


%-

lme_juice_fitted = fitted(lme_juice)
%ci = coefCI(lme_juice)
figure;
for ar = 1 : length(area2test)
    take = ismember(modeldata.area,area2test{ar})
    plot(tbl_juice.pref(take),lme_juice_fitted(take),'o','Color',colorsArea(ar,:));
    hold on
end

lme_pb12_fitted = fitted(lme_pb12)
%ci = coefCI(lme_juice)
figure;
for ar = 1 : length(area2test)
    take = ismember(modeldata.area,area2test{ar})
    plot(tbl_pb12.pref(take),lme_pb12_fitted(take),'o','Color',colorsArea(ar,:));
    hold on
end



%% FIG S4 - individual monkeys!

grp = {'modeldata.perf_ju' , 'abs(modeldata.pref-0.5)' ; ...
    'modeldata.perf_pb' , 'abs(modeldata.pref-0.5)' ;...
    'modeldata.perf_unpb' , 'abs(modeldata.pref-0.5)' ;...
    'modeldata.perf_si' , 'abs(modeldata.pref-0.5)' ;...
    'modeldata.perf_pb2-modeldata.perf_pb1' , 'modeldata.pref' }
yscs = [0.4 0.9 ; 0.15 0.55 ;0.15 0.55 ; 0.4 0.9 ; -0.15 0.15];

colorsPr = cbrewer('seq', 'Greys', 100);
figure;
sub = 0;
for mk = 1 : 2
    for m = 1 : size(grp,1)
        y_all = eval([grp{m,1}]) ;
        x_all = eval([grp{m,2}]) ;

        ysc = [round(min(y_all),2)-.05 round(max(y_all),2)+.05];
        ysc = yscs(m,:);
        if strcmp(grp{m,2},'abs(modeldata.pref-0.5)')
            xsc = 0.5;
        else
            xsc = 1;
        end
        sub = sub + 1;
        subplot(2,size(grp,1),sub)

        for ar = 1 : length(area2test)
            %   take = ismember(modeldata.area,area2test(ar)) & ~isnan(x_all) & ~isnan(y_all);
            if mk==1
                take = ismember(modeldata.area,area2test(ar)) & ismember(modeldata.mk,'M') & ~isnan(x_all) & ~isnan(y_all);
            else
                take = ismember(modeldata.area,area2test(ar)) & ismember(modeldata.mk,'X') & ~isnan(x_all) & ~isnan(y_all);

            end
            x = x_all(take);
            y = y_all(take);
            xsc = max(x);
            [r,p] = corrcoef(x,y)
            allRP(ar,:) = [r(2,1) , p(2,1)];

            xgrid=linspace(min(x),xsc,100);
            ygrid=linspace(ysc(1),ysc(1),100);
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

            colormap(colorsPr)
            hold on
            line([0.5 0.5],[ysc(1) ysc(2)],'Color',[.6 .6 .6])
            line([0 1],[0 0],'Color',[.6 .6 .6])
            %   plot(x,y,'.','Color',colorsArea_sub(ar,:))
            text(.1,ysc(1)+0.02+(ar/70),['\it r(' num2str(length(x)-2) ')=' num2str(round(allRP(ar,1),3)) ', p=' num2str(round_pval(allRP(ar,2)))],...
                'Color',colorsArea(ar,:))
            %  text(xsc/2+.1,ysc(1)+0.02+(ar/70),['p=' num2str(round_pval(allRP(ar,2)))],'Color',colorsArea(ar,:))

            cf = fit(x,y,'poly1'); % fit
            % p_handle = plot(cf,'k'); % Plot fit
            p_handle = plot(xgrid,cf.p1*xgrid + cf.p2)

            if  allRP(ar,2)<0.05
                set(p_handle,'Color',colorsArea(ar,:),'LineWidth',3);
            else
                set(p_handle,'Color',colorsArea(ar,:),'LineWidth',1,'LineStyle','-');
            end
            legend off
            xlabel(grp{m,2});ylabel(grp{m,1})
            % title(area2test{ar})
            if strcmp(grp{m,2},'abs(modeldata.pref-0.5)')
                xlim([0 .5]);
            else
                xlim([0 1]);
            end
            ylim([ysc(1) ysc(2)])
            box on

        end
    end
end

