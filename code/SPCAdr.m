function[var_set,rmse_SPCA, Fail, Theta, R2_OSS_SPCA,Nr,spca_out] = SPCAdr(troll,Nfacmax,yfor,X,yreg,N)
% Rolling Sample Forecasting
% Supervised PCA
% Nfacmax          % Number of factors max
% troll         % Length of the rolling window 
% theta0        % start point of interval for theta to consider under cross-calidation
% theta1        end point of interval for theta to consider under cross-calidation
% theta0        % start point of interval for theta to consider under cross-calidation
% DELTA         % size of step through the interval

%%
theta0 = 0.01;    
DELTA = 0.001;      
theta1 = 0.02;    

xst = normalize(X(1:troll,:),1); 
yst = normalize(yfor(1:troll));        

 
%% Supervision: measure the correlation between y and the
    % explanatory variables
    % Compute standard univariate regression coefficients for each variable
    % Calculate the score (either correlation or p-value)
  g = zeros(1,N);
  for j = 1:N
      g(1,j) = (xst(:,j)'*xst(:,j))\(xst(:,j)'*yst);   
      if ismissing(g(1,j))== 1
        g(1,j)=0; 
      end
  end
%% Cross-validation for determining the supervised PCA threshold
rmse = [];
thetaN = 0;
spca_out = NaN(length(theta0:(DELTA):theta1),size(yfor,1),Nfacmax);
for theta = theta0:(DELTA):theta1
    disp("Theta loop " +theta)
    thetaN = thetaN + 1;
    %pick the variables to keep  
    IND = zeros(size(xst));
    for j = 1:N  
        if abs(g(1,j)) > theta
            IND(:,j)=1;
        end
    end
                            
    if theta == theta0
       var_set = find(sum(IND.*xst));
    elseif size(find(sum(IND.*xst)),2)==size(var_set,2)
       continue  
    else
       var_set = find(sum(IND.*xst)); % the set of explanatory variables for PCA
    end
    if size(var_set,2)>Nfacmax
        [B,var_set_reg] = maxk(abs(IND(1,:).*g),Nfacmax); %picks the Nfacmax number of variables, which are most strongly correlated with y
    end
                     
    %% Carry out Forecasting over Rolling Samples %
    % Matrices for saving forecasts % 
    yfcst_spc = NaN(size(yfor,1),Nfacmax);  %OLS with supervised PCA factors
    yfcst_reg = NaN(size(yfor,1),1);
    yfcst_ar = NaN(size(yfor,1),1);
    ymean = NaN(size(yfor,1),1);   
    
    %%Starting the rolling window loop for estmation
    Nr=NaN(size(X,1),1); %number of variables included in the regression (in case of multicollinearity can be smaller than Nf)
    for t = 1:1:size(X,1)-troll
        % Check if there are any time gaps in the considered window
        ylag_t =yreg(t:t+troll,2);
        gap = find(isnan(ylag_t(1:end,:)), 1);
        if isempty(gap) == 0 %roll the window one unit further in case 
           continue            
        end
        
        % Take the rolling window from the full sample
        xst_t = X(t:t+troll,:);
        yfor_t = yfor(t:t+troll);
        yreg_t = yreg(t:t+troll,1:1);
    
        % Lag Series to Eliminate Leads of y in yfor %
        yfor_t_lag = yfor_t(1:troll);
        yreg_t_lag = yreg_t(1:troll,:);
        ylag_t_lag=ylag_t(end-troll:end-1,:);
                       
       
        %% % Estimate a historical mean
        ymean(troll+t) = mean(rmmissing(ylag_t_lag));
    
        %E estimating a model with a constant only
        y = yfor_t_lag;
        x = yreg_t_lag;
        beta_ar = x\y;       
        yfcst_ar(troll+t) = yreg_t(end,1)*beta_ar;
         
        clear x;       
           
        %% %pick the subset of explanatory variables for PCA
          xst_t_st=rmmissing(xst_t(:,var_set),2); 
    
         % Standardise the explanatory variables prior to running PCA
           xst_t_st = normalize(xst_t_st);
           xst_t_st = rmmissing(xst_t_st,2);
    
        %% %Estimate a straightforward regression on the chosen subset of
        % variables if too many selected variables, limit to 
        if size(var_set,2)>Nfacmax
             xst_reg_t = rmmissing(normalize(xst_t(:,var_set_reg)),2);
             NR=size(xst_reg_t,2); %number of variables included in the regression (in case of multicollinearity can be smaller than Nf)
             J=1;
             while J <=NR
                if var(xst_reg_t(:,J))==0
                   xst_reg_t(:,J)=[];
                   NR=size(xst_reg_t,2);
                end
                J=J+1;
             end
             Nr(t,1)=NR;
             x =[ones(troll,1) xst_reg_t(1:troll,:)];
             b_reg = inv(x'*x)*x'*y;
             
             yfcst_reg(troll+t,1) = [1 xst_reg_t(end,:)]*b_reg; 
        else
           %Estimate a straightforward regression on the chosen subset of
           %variables if reasonable number of variables in set         
              x =[ones(troll,1) xst_t_st(1:troll,:)];
              b_reg = (x'*x)\x'*y;
    
             yfcst_reg(troll+t,1) =[1 xst_t_st(end,:)]*b_reg;  
        end   
        clear x;
         %% %Obtain the factor estimates
           Nfmax = min(size(xst_t_st,2),Nfacmax); % limit the number of future components in case the subset is too small
           [~,Fhat] = pca(xst_t_st,'NumComponents',Nfmax,'Centered',false);
           % Estimate OLS with PC factors, with varying number of components               
           y = yfor_t_lag; 
           for jj = 1:Nfmax
               Fhat_est = Fhat(:,1:jj);
               Fhat_x =[yreg_t_lag Fhat_est(1:troll,:)];
               beta_ols = (Fhat_x'*y)/size(Fhat_x,1);
               yfcst_spc(troll+t,jj) = [yreg_t(end,1) Fhat_est(end,:)]*beta_ols; 
           end
    end
    
    yfcst_all = [ymean, yfcst_reg, yfcst_spc(:,1:Nfmax)];
    spca_out(thetaN,:,1:size(yfcst_spc,2)) = yfcst_spc(:,1:Nfacmax);
    % Construct residuals %
    E = rmmissing(yfor - yfcst_all,1);
        Fail = 0;
        if size(rmmissing(E(:,1)),1) < size(rmmissing(E(:,2)),1)
           Fail = 1; 
        end
    rmse_e=zeros(1,size(E,2));
    for k = 1:size(E,2)
        e =rmmissing(E(:,k));       
        rmse_e(1,k) = sqrt(mean(e.^2,1));    
    end    
    
    if isempty(rmse)==1
        rmse = rmse_e; %RMSE*10^3 measure of forecast combinations 
    else
        out_rm = rmse_e;      
        L=max(numel(out_rm),numel(rmse(1,:)));
        out_rm(:,end+1:L)=nan;
        rmse(:,(end+1:L))=nan;
        rmse = [rmse;out_rm];
    end

    
end % Theta loop ends
[row,col] = find((rmse(:,3:end)*10^3)== min(rmse(:,3)*10^3,[],'all'));
rmse_SPCA = rmmissing([rmse(row,:)]);
spca_out = squeeze(spca_out(row,:,:));
Theta = [theta0:(DELTA):theta1]';
Theta = Theta(row);
Nr=mean(rmmissing(Nr));
R2_OSS_SPCA = 100*(1-(rmse_SPCA(1,2:end).^2/rmse_SPCA(1,1).^2)); %Campbell-Thompson out-of-sample R^2
end 