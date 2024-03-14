function [perf,Post_all,out] = pca_lda_kfold(XX,YY,param)

ncd = length(unique(YY{1}));
Post_all = zeros(length(XX),ncd,ncd);
for t = 1 : length(XX)
    group=YY{t};
    data_sub = XX{t};
    
    [coeff,score,~,~,explained,~] = pca(data_sub,'NumComponents',param.nComp);
    data_sub = score;

    cv = cvpartition(group,'KFold',10); % 0.20'
    numFolds = cv.NumTestSets;
    for i = 1:numFolds
        sampleInds = cv.test(i);
        trainInds = cv.training(i);

        trainingData = data_sub(trainInds,:);
        sampleData = data_sub(sampleInds,:);
        
        %- try classification.. failed when too many units have non-zeros FR
        class = [];
        try [class,err,posterior,logp,coeff] = classify(sampleData,trainingData,group(trainInds), 'diaglinear');
        end

        if ~isempty(class)
            testi = group(sampleInds);
            cd=unique(group);
            for tr = 1 : length(class)
                Post_all(t,cd==class(tr),cd==testi(tr)) = Post_all(t,cd==class(tr),cd==testi(tr))+1;
            end
            mean(abs(class-group(sampleInds)));
            perf(t,i) = mean(abs(class-group(sampleInds))==0);
        else
            perf(t,i) = NaN;
        end

        
    end
end
perf = nanmean(perf,2);

for t = 1 : length(XX)
    post_dumm = squeeze(Post_all(t,:,:));
    Cds = 1:size(post_dumm,1);
    for cl = 1 : length(Cds)
        TP = post_dumm(cl,cl);
        FN = sum(post_dumm(cl,Cds~=cl));
        FP = sum(post_dumm(Cds~=cl,cl));
        TN = sum(sum(post_dumm(Cds~=cl,Cds~=cl)));

        out.recall(t,cl) = TP/(TP+FN);
        out.precision(t,cl) = TP/(TP+FP);
        out.specificity(t,cl) = TN/(TN+FP);
    end
end
out.f_meas = (2.*out.recall.*out.precision) ./ (out.recall+out.precision)   ;


