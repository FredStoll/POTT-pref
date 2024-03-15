%% POTT Pref - Model comparison
%-
%-
%- Author: Fred M. Stoll, Icahn School of Medicine at Mount Sinai, NY
%- Date: 2022.12

clear

%- locate the files
path2go = '/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-BehavPrefChange/data-final/';

skip = true; %- re-run the models or just run on saved models!
% figure;
if ~skip
    mk = 'Morbier'; % 'Morbier'

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
    x=0;
    for sess = 1 : length(list) % for every sessions
        filename = [path2go list(sess).name];
        disp(sess)
        [ALLtemp,param] = behav_model_revision(filename,false,showfig); %- no permutations
        if ~isempty(ALLtemp)
            x = x + 1;
            ALL(x,1) = ALLtemp;
        end
        if perm
            for p = 1 : nPerm
                [ALLperm(sess,p),~] = behav_model_revision(filename,perm,showfig);
            end
        end

    end
    if perm
        save([path2go mk '_behav_revision.mat'],'ALL','ALLperm','param','-v7.3')
    else
        save([path2go mk '_behav_revision.mat'],'ALL','param','-v7.3')
    end
end

%% posthoc - compare proportion of sessions better explained by each model (juice agnostic or not)

clear
%- locate the files
path2go = '/Users/fred/Dropbox/Rudebeck Lab/ANA-POTT-BehavPrefChange/data-final/';

M1 = load([path2go 'Morbier_behav_revision.mat'],'ALL','param')
M2 = load([path2go 'Mimic_behav_revision.mat'],'ALL','param')

%- remove sessions before threshold for Morbier, if needed
subset = {'052318' '010122'} ; %- for anything after last change
if ~isempty(subset)
    lim = datenum(subset,'mmddyy');
    for ss = 1 : length(M1.ALL)
        if ~isempty(M1.ALL(ss).name)
        days(ss) = datenum(M1.ALL(ss).name(end-16:end-11),'mmddyy');
        else 
            days(ss)=0
        end
    end
    takeme = (days>=lim(1) & days<lim(2)) ;
else
    takeme = true(1,length(M1.ALL));
end
M1.ALL = M1.ALL(takeme);

%- number trials per sessions (all set, not only converging ones)
AIC=[];
pref=[];
mk=[];
r2=[];
session_name={};
for m = 1 : 2
    eval(['ALL = M' num2str(m) '.ALL;']);
    for s = 1 : length(ALL)
        if ALL(s).converge_mdl1 & ALL(s).converge_mdl2 
            AIC = [AIC ; ALL(s).mdl1.ModelCriterion.AIC ALL(s).mdl2.ModelCriterion.AIC];
            r2 = [r2 ; ALL(s).mdl1.Rsquared.Adjusted ALL(s).mdl2.Rsquared.Adjusted];
            pref = [pref ; ALL(s).pref];
            mk = [mk ; m];
            session_name = [session_name ; ALL(s).name(end-17:end-10)];
        end
    end
end

prop_mdl = [sum(sign(AIC(:,1)-AIC(:,2))==-1) length(AIC) ;...
            sum(sign(AIC(:,1)-AIC(:,2))==-1 & mk==1) length(AIC(mk==1)) ;...   
            sum(sign(AIC(:,1)-AIC(:,2))==-1 & mk==2) length(AIC(mk==2))]  ;

prop_mdl(:,1)./prop_mdl(:,2)

for m = 1 : 3
    [tbl,chi2stat(m),pval(m)] = chi2_fms(prop_mdl(m,1),prop_mdl(m,2),prop_mdl(m,2)-prop_mdl(m,1),prop_mdl(m,2))
end

figure;histogram(abs(pref(sign(AIC(:,1)-AIC(:,2))==-1)-0.5),[0:0.025:0.5]);
hold on
histogram(abs(pref(sign(AIC(:,1)-AIC(:,2))==1)-0.5),[0:0.025:0.5]);
xlabel('Absolute Preference');
ylabel('Number of sessions')

keep_sessions = session_name(sign(AIC(:,1)-AIC(:,2))==1);

save([path2go 'POTT_Behav_subset.mat'],'keep_sessions')

prctile(abs(pref(sign(AIC(:,1)-AIC(:,2))==-1)-0.5),[50 25 75])
prctile(abs(pref(sign(AIC(:,1)-AIC(:,2))==1)-0.5),[50 25 75])

kwtest = [abs(pref(sign(AIC(:,1)-AIC(:,2))==-1)-0.5), ones(size(abs(pref(sign(AIC(:,1)-AIC(:,2))==-1)-0.5))) ; ...
          abs(pref(sign(AIC(:,1)-AIC(:,2))==1)-0.5), 2*ones(size(abs(pref(sign(AIC(:,1)-AIC(:,2))==1)-0.5)))];

kruskalwallis(kwtest(:,1),kwtest(:,2))
