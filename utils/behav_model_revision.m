
function [ALL,param] = behav_model_revision(filename,perm,showfig)

%- load behav variable
load(filename)

%- some utils
take = @(arg) strcmp(TrialType_header,arg);
norm = @(data) -1+((data-min(data))*2)/(max(data)-min(data)) ;
when = @(arg,cd) TrialType(strcmp(TrialType_header,arg),:)==cd ;

%- permutations
if perm
    
    pav = when('task',1); %- leave the PAV trial alone
    ntr4perm = sum(~pav);
    tr_order = randperm(ntr4perm); %- perm the others (ins + dyn)
    
    param2perm = {'TrialType' 't_FPon' 't_FPfix_first' 't_respfix' 't_resp'};
    for pp = 1 : length(param2perm)
        eval(['dumm = ' param2perm{pp} '(:,~pav);']);
        eval([param2perm{pp} '=[' param2perm{pp} '(:,pav),dumm(:,tr_order)];']);
    end
    
    %- reboot the function to send the new TrialType matrix!
    when = @(arg,cd) TrialType(strcmp(TrialType_header,arg),:)==cd ;
end


%- trials to consider
diff_juice = (when('I_juiceL',1) & when('I_juiceR',2)) | (when('I_juiceL',2) & when('I_juiceR',1)); %- take only the different juice trials
restr = when('task',2) & when('brks',0) & diff_juice; %- take only INS trials (task == 2), completed (brks == 0) and diff juice

%- variable to predict: ChooseJ1 - 1 for J1 / 0 for J2
choice = abs(TrialType(take('I_chosenjuice'),restr)-2);

if length(choice) > 100  %- trash sessions with less than 100 trials diff juice
    
    %% Step 1 - GLM on Choices in different juice trials
    %- extract Proba of both options
    probaJ1( when('I_juiceL',1) & restr ) = TrialType(take('I_probaL'),when('I_juiceL',1) & restr);
    probaJ1( when('I_juiceR',1) & restr ) = TrialType(take('I_probaR'),when('I_juiceR',1) & restr);
    probaJ1 = probaJ1(restr);
    probaJ2( when('I_juiceL',2) & restr ) = TrialType(take('I_probaL'),when('I_juiceL',2) & restr);
    probaJ2( when('I_juiceR',2) & restr ) = TrialType(take('I_probaR'),when('I_juiceR',2) & restr);
    probaJ2 = probaJ2(restr);
    
    pbJ1 = norm(probaJ1/100); %- normalize proba
    pbJ2 = norm(probaJ2/100);
    
    %- extract Reaction time
    rt = t_respfix-t_resp(1,:);
    rt = norm(log(rt(restr))); %- log then min-max normalization (-1 to 1)
    rt(~isfinite(rt)) = NaN;
    
    %- extract Initiation time
    ft = t_FPfix_first-t_FPon(1,:);
    ft = norm(log(ft(restr))); %- log then min-max normalization (-1 to 1)
    ft(~isfinite(ft)) = NaN;
    
    %- extract prevJuice (0 if chose J1 before / 1 if chose J2 before) - this is not side related
    %- prevJuice == which juice he 'could have get' on the last completed trial (not necessarily
    %- rewarded).
    %- WARNING!!! FOR NOW, INCLUDE Dynamic task trials -- If trial before is Dynamic task, which juice he chose there (Morbier) or was offered (Mimic)
    
    currjuice = NaN(1,length(TrialType));
    currjuice(when('task',1)) = TrialType(take('P_juice'),when('task',1))-1;
    currjuice(when('task',2)) = TrialType(take('I_chosenjuice'),when('task',2))-1;
    currjuice(when('task',5)) = TrialType(take('D_chosenjuice'),when('task',5))-1;
    currjuice(when('task',3)) = TrialType(take('D_chosenjuice'),when('task',3))-1;
    
    %- put the last trials' juice when brks at current (or multiple breaks in a row, take the last one)
    noJuice = find(isnan(currjuice));
    for tr = 1 : length(noJuice)
        if noJuice(tr)>1
            currjuice(noJuice(tr)) = currjuice(noJuice(tr)-1);
        end
    end
    
    prevJuice = [NaN currjuice(1:end-1)];
    prevJuice = prevJuice(restr);
    prevJuice(prevJuice==0)=-1; %- a way to standardize it (-1 vs 1)
    
    %- extract the last received juice
    received = when('rew',1);
    currjuice(~received)=NaN; %- remove all the juice info from unrewarded trials
    
    %- put the last trials' RECEIVED juice when no rew or brks (or multiple breaks in a row, take the last one)
    %- no information on how far was this last reward
    noRwd = find(isnan(currjuice));
    for tr = 1 : length(noRwd)
        if noRwd(tr)>1
            currjuice(noRwd(tr)) = currjuice(noRwd(tr)-1);
        end
    end
    
    prevJuiceReceived = [NaN currjuice(1:end-1)];
    prevJuiceReceived = prevJuiceReceived(restr);
    prevJuiceReceived(prevJuiceReceived==0)=-1; %- a way to standardize it (-1 vs 1)
    
    modelspec1 = 'choice ~ 1 + probDiff + prevJuice' ; %- model 1 choice only based on comparison of proba
    modelspec2 = 'choice ~ 1 + probJ1 + probJ2 + prevJuice' ;   %- model 2 choice based on both proba and juice

    T = table(choice',pbJ1',pbJ2',(pbJ1-pbJ2)',ft',rt',prevJuice',prevJuiceReceived','VariableNames',{'choice','probJ1','probJ2','probDiff','ft','rt','prevJuice','prevJuiceReceived'});
    T.choice=categorical(T.choice);
    T.prevJuice=categorical(T.prevJuice);
    T.prevJuiceReceived=categorical(T.prevJuiceReceived);
    
    %- 1st model
    lastwarn('', '');
    mdl1 = fitglm(T,modelspec1,'Distribution','binomial','Link','logit');
    %- check if glm converged
    [~, warnId] = lastwarn();
    if strcmp(warnId,'stats:glmfit:IterationLimit')
        converge_mdl1=0;
    else
        converge_mdl1=1;
    end
    
    lastwarn('', '');
    mdl2 = fitglm(T,modelspec2,'Distribution','binomial','Link','logit');
    %- check if glm converged
    [~, warnId] = lastwarn();
    if strcmp(warnId,'stats:glmfit:IterationLimit')
        converge_mdl2=0;
    else
        converge_mdl2=1;
    end

    disp([mdl1.ModelCriterion.BIC mdl2.ModelCriterion.BIC])
    disp([mdl1.Rsquared.Adjusted mdl2.Rsquared.Adjusted])

    %- predict choices (for that, prevJuice is set as 0 (J1) and RT is median(rt)
    clear newf newc
    probas = norm([10:1:90]);
    thr = NaN(length(probas),1);
    
    for i = 1 : length(probas)
        
        tab = table((probas(i)*ones(length(probas),1)),probas',nanmedian(ft)*ones(length(probas),1),nanmedian(rt)*ones(length(probas),1),-ones(length(probas),1),-ones(length(probas),1),'VariableNames',{'probJ1','probJ2','ft','rt','prevJuice','prevJuiceReceived'});
        tab.prevJuice = categorical(tab.prevJuice);
        tab.prevJuiceReceived = categorical(tab.prevJuiceReceived);
        
        [newf(i,:) , newc] = predict(mdl2,tab); %- last column is the sessiooon
        
        if ~isempty(find(newf(i,:)>=0.5,1,'first'))
            thr(i,1) = find(newf(i,:)>=0.5,1,'first');
        end
    end
    
    %- extract the preference for J1 vs J2
    %- sensitivity doesn't take into account the interaction factor... so can use the nb of
    %- estimated probaChoiceJ1>0.5.. if pref>0.5, prefer the J1
    sensitivity = mdl2.Coefficients{'probJ1','Estimate'}/mdl2.Coefficients{'probJ2','Estimate'};
    pref = sum(sum(newf<0.5))/(length(newf(:))-sum(sum(newf==0.5)));

    %% Extract param of interest for each session!
    ALL.name = filename;
    ALL.mdl1 = mdl1;
    ALL.converge_mdl1 = converge_mdl1;
    ALL.mdl2 = mdl2;
    ALL.converge_mdl2 = converge_mdl2;
    ALL.pref = pref;

    if ~perm %- keep more param when not perm
        ALL.TrialType = single(TrialType);
        ALL.trial2take = restr;
        ALL.T = T;
        ALL.sensitivity = sensitivity;
        ALL.rt = single(rt);
        ALL.ft = single(ft);
   end
    
    param.modelspec = {modelspec1 modelspec2};

else
    ALL = [];
    param = [];
end

