function ci = get_bootstrap_cis(alpha, stat, bstat, bootfun, data, type, weights)
        
% convenience function modified from MATLAB's bootci.m 
% Allows for calculation of different kinds of CIs at different alpha levels with pre-generated bootstrap statistics (i.e., bootci makes you rerun the bootstrapping every time)
% Joseph Griffis 2024

bootstrpOptions = statset('bootstrp');

switch type
    case 'bcp'
        % Given bootstrapped samples, get BCa CIs at desired alpha
    
        z_0 = fz0(bstat,stat);
        z_alpha = norminv(alpha/2); % normal confidence point
         
        % transform z0 back using the invECDF[normCDF(2z0-za)] and
        % invECDF[normCDF(2z0+za)] 
        pct1 = 100*normcdf(2*z_0-z_alpha); 
        pct2 = 100*normcdf(2*z_0+z_alpha);
        
        % inverse of ECDF
        m = size(bstat,2);
        lower = zeros(1,m);
        upper = zeros(1,m);
        for i=1:m
            lower(i) = prctile(bstat(:,i),pct2(i),1);
            upper(i) = prctile(bstat(:,i),pct1(i),1);
        end
        ci = [lower;upper];

    case 'bca'
        % same as bootcper, this is the bias correction
        z_0 = fz0(bstat,stat);
        
        % apply jackknife
        try
            jstat = jackknife(bootfun,data{:},'Options',bootstrpOptions);
        catch ME
            m = message('stats:bootci:JackknifeFailed',func2str(bootfun));
            throw(addCause(MException(m.Identifier,'%s',getString(m)),ME));
        end
        N = size(jstat,1);
        if isempty(weights)
            weights = repmat(1/N,N,1);
        else
            weights = weights(:);
            weights = weights/sum(weights);
        end
        
        % acceleration finding, see DiCiccio and Efron (1996)
        mjstat = sum(jstat.*weights,1); % mean along 1st dim.
        score = mjstat - jstat; % score function at stat; ignore (N-1) factor because it cancels out in the skew
        iszer = all(score==0,1);
        skew = sum((score.^3).*weights,1) ./ ...
            (sum((score.^2).*weights,1).^1.5) /sqrt(N); % skewness of the score function
        skew(iszer) = 0;
        acc = skew/6;  % acceleration
        
        % transform back with bias corrected and acceleration
        z_alpha1 = norminv(alpha/2);
        z_alpha2 = -z_alpha1;
        pct1 = 100*normcdf(z_0 +(z_0+z_alpha1)./(1-acc.*(z_0+z_alpha1)));
        pct1(z_0==Inf) = 100;
        pct1(z_0==-Inf) = 0;
        pct2 = 100*normcdf(z_0 +(z_0+z_alpha2)./(1-acc.*(z_0+z_alpha2)));
        pct2(z_0==Inf) = 100;
        pct2(z_0==-Inf) = 0;
        
        % inverse of ECDF
        m = numel(stat);
        lower = zeros(1,m);
        upper = zeros(1,m);
        for i=1:m
            lower(i) = prctile(bstat(:,i),pct2(i),1);
            upper(i) = prctile(bstat(:,i),pct1(i),1);
        end
        
        % return
        ci = sort([lower;upper],1);

    case 'normal'
        se = std(bstat,0,1);   % standard deviation estimate
        bias = mean(bstat - stat,1);
        za = norminv(alpha/2);   % normal confidence point
        lower = stat - bias + se*za; % lower bound
        upper = stat - bias - se*za;  % upper bound
        
        % return
        ci = [lower;upper];        
end

end

function z0=fz0(bstat,stat)
% Compute bias-correction constant z0
z0 = norminv(mean(bstat < stat,1) + mean(bstat == stat,1)/2);
end   