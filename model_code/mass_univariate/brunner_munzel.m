function [testStat, p_value] = brunner_munzel(sample1, sample2)
    
    %some initialization
    numSamples1 = length(sample1);
    numSamples2 = length(sample2);
    N = numSamples1 + numSamples2;
    rankSum1 = 0;
    rankSum2 = 0;
   
    %get mean ranksum of sample1
    for k=1:numSamples1 
        tmpRank = 0.5;
        for j=1:2
            if j==1
                numSamples = numSamples1;
                for l=1:numSamples
                    if sample1(k)-sample1(l) < 0
                        c = 0;
                    elseif sample1(k)-sample1(l) == 0
                        c = 0.5;
                    else 
                        c = 1;
                    end
                    tmpRank = tmpRank + c;
                end   
            else
                numSamples = numSamples2;
                for l=1:numSamples
                    if sample1(k)-sample2(l) < 0
                        c = 0;
                    elseif sample1(k)-sample2(l) == 0
                        c = 0.5;
                    else 
                        c = 1;
                    end
                    tmpRank = tmpRank + c;
                end   
            end
        end
        rankSum1 = rankSum1 + (1/numSamples1)*tmpRank ;
    end

    %get mean ranksum 2
    for k=1:numSamples2 
        tmpRank = 0.5;
        for j=1:2
            if j==1
                numSamples = numSamples1;
                for l=1:numSamples
                    if sample2(k)-sample1(l) < 0
                        c = 0;
                    elseif sample2(k)-sample1(l) == 0
                        c = 0.5;
                    else 
                        c = 1;
                    end
                    
                    tmpRank = tmpRank + c;
                end   
            else
                numSamples = numSamples2;
                for l=1:numSamples
                    if sample2(k)-sample2(l) < 0
                        c = 0;
                    elseif sample2(k)-sample2(l) == 0
                        c = 0.5;
                    else 
                        c = 1;
                    end
                    
                    tmpRank = tmpRank + c;
                end   
            end
        end
        rankSum2 = rankSum2 + (1/numSamples2)*tmpRank ;
    end

    %calculate variance S1_squared
    S1sq = 0;
    for k=1:numSamples1 
        tmpRank = 0.5;
        internRank1 = 0.5;
        for j=1:2
            if j==1
                numSamples = numSamples1;
                for l=1:numSamples
                    if sample1(k)-sample1(l) < 0
                        c = 0;
                    elseif sample1(k)-sample1(l) == 0
                        c = 0.5;
                    else 
                        c = 1;
                    end
                    tmpRank = tmpRank + c;
                    internRank1 = internRank1 + c;
                end   
            else
                numSamples = numSamples2;
                for l=1:numSamples
                    if sample1(k)-sample2(l) < 0
                        c = 0;
                    elseif sample1(k)-sample2(l) == 0
                        c = 0.5;
                    else 
                        c = 1;
                    end
                    
                    tmpRank = tmpRank + c;
                end   
            end
        end
        S1sq = S1sq + (1/(numSamples1-1))*(tmpRank - internRank1 - rankSum1 + (numSamples1 +1)/2)^2;
    end
    
    %calculate variance S_squared
    S2sq = 0;
    for k=1:numSamples2 
        tmpRank = 0.5;
        internRank2 = 0.5;
        for j=1:2
            if j==1
                numSamples = numSamples1;
                for l=1:numSamples
                    if sample2(k)-sample1(l) < 0
                        c = 0;
                    elseif sample2(k)-sample1(l) == 0
                        c = 0.5;
                    else 
                        c = 1;
                    end
                    tmpRank = tmpRank + c;
                    
                end   
            else
                numSamples = numSamples2;
                for l=1:numSamples
                    if sample2(k)-sample2(l) < 0
                        c = 0;
                    elseif sample2(k)-sample2(l) == 0
                        c = 0.5;
                    else 
                        c = 1;
                    end
                    
                    tmpRank = tmpRank + c;
                    internRank2 = internRank2 + c;
                end   
            end
        end
        S2sq = S2sq + (1/(numSamples2-1))*(tmpRank - internRank2 - rankSum2 + (numSamples2 +1)/2)^2;
    end
    
 
    %final variance
    var = N*S1sq/(N-numSamples1) + N*S2sq/(N-numSamples2);
    
    %calcuate the value of the test statistic
    testStat = ((rankSum2-rankSum1)/(sqrt(var)))*sqrt(numSamples1*numSamples2/N);
    
    %calculate p-value
    p_value =  2*min(1-normcdf(testStat), normcdf(testStat));
   
end
