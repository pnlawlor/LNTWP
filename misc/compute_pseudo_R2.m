function R2 = compute_pseudo_R2(Y, Yhat, Ymean, distr, Yhat_null)

if(nargin < 4)
  Modelr = sum(Y.*log(eps+Yhat) - Yhat);
  Intercr = sum(Y.*log(eps+Ymean) - Ymean);
  Sat_r = sum(Y.*log(eps+Y) - Y);


  R2 = (1-(Sat_r-Modelr)./(Sat_r-Intercr))';
end

if(nargin >= 4 && strcmp(distr, 'logit'))
  % Log likelihood of model under consideration
  L1 = 2*length(Y)*sum(Y.*log((Yhat==0)+Yhat)./mean(Yhat) + (1-Y).*log((Yhat==1)+1-Yhat)./(1-mean(Yhat)));
  
  % Log likelihood of homogeneous model
  %b0_hat_null = mean(log((Y==0)+Y) - log((Y==1)+1-Y));
  
  %b0_hat_null = glmfit(ones(length(Y),1), Y, 'binomial', 'link', 'logit', 'constant', 'off');
  %Yhat_null = 1/(1+exp(-b0_hat_null));
  
  %Yhat_null = mean(Yhat);
  
  L0 = 2*length(Y)*sum(Y.*log((Yhat_null==0)+Yhat_null)./mean(Yhat) + (1-Y).*log((Yhat_null==1)+1-Yhat_null)./(1-mean(Yhat)));
  
  % Log likelihood of saturated model
  Lsat = 2*length(Y)*sum(Y.*log((Y==0)+Y)./mean(Y) + (1-Y).*log((Y==1)+1-Y)./(1-mean(Y)));
  
  % Note that saturated log likelihood is 0 for binomial distribution
  %R2 = (1-(Lsat - L1)./(Lsat-L0));
  R2 = 1 - L1/L0;
end

if(nargin >= 4 && strcmp(distr, 'gamma'))
    k = (mean(Y).^2./var(Y));
    theta = var(Y)./mean(Y);
    
    % Log likelihood of model under consideration
    L1 = (k-1)*sum(log(eps+Yhat)) - sum(Yhat)/theta;
  
    % Log likelihood of homogeneous model
    L0 = (k-1)*length(Y)*log(eps+Ymean) - length(Y)*Ymean/theta;
    
    % Log likelihood of saturated model
    Lsat = (k-1)*sum(log(eps+Y)) - sum(Y)/theta;
    R2 = (1-(Lsat - L1)./(Lsat-L0));
end

if(nargin >= 4 && (strcmp(distr, 'normal') || strcmp(distr, 'gaussian')))
    R2 = 1 - sum((Y-Yhat).^2)/sum((Y-Ymean).^2);
end