
function [] = main_function_replication(xdata,md_Y,cases,troll,frequency,savepath)
    % Estimation
    est_pca = 1;
    est_spca = 1;
    est_lasso = 1;
    est_ada_lasso = 1;
    est_rf = 1;
    Nfacmax = 5;
    when_end = 201037; % 201037 is whole sample
    rep = 1; % Sets seed value so CV produces reproducable results in Lasso and RF
    
    % Re-estimation frequency
    if contains('minutely',frequency) % Hourly frequency
        re_estimate_lasso = floor(troll/5); % How  many obs use before re-estimation (to save time)
        re_estimate_rf = floor(troll/5); % How  many obs use before re-estimation (to save time)
    else
        re_estimate_lasso = floor(troll/20); % How  many obs use before re-estimation (to save time)
        re_estimate_rf = floor(troll/20); % How  many obs use before re-estimation (to save time)
    end
    xdata = xdata(1:when_end,:);
    md_Y = md_Y(1:when_end,:);
    X = table2array(xdata);
    
    [~,N] = size(X);
    %% Handle missing values
    for j = 1:N
        X(isnan(X(:,j)),j) = mean(X(:,j),'omitnan');
    end
    X = timetable(xdata.Time,X);
    % Put X and y in same timetable
    Xy = innerjoin(X,md_Y);
    
    % Change frequency
    if contains('hourly',frequency)
        Xy = retime(Xy,frequency,'sum'); 
    end

    % remove potential zeros added by retiming
    allZeros = all(Xy.Variables == 0,2);
    Xy(allZeros,:) = [];
  
    X = table2array(Xy(:,1:end-5));
    Y = table2array(Xy(:,end-4:end));
    [T,N] = size(X);
    Nf = min(N,Nfacmax);
    %% From this point forward depends on currency 
    all_currencies = {'USDCAD','AUDUSD','GBPUSD','USDJPY','EURUSD'};
    clear md* xdata
    VAR_OUT = [];
    RMSE_SPCA = nan(5,12);
    R2OS_SPCA = nan(5,10);
    % empty variables so save works even if they are not estimated
    save_coeffs_lasso = [];
    save_coeffs_adalasso = [];
    for currency = 1:5
        curr = all_currencies{currency};
        disp("Working on currency " +curr)
        y = Y(:,currency);
        yfor=[y(2:end,1);nan]; %generating the lead of the explanatory variable to be used in the predictive regression
        
        if cases == 4 || cases == 10 % For this case X depends on currency
            Xtmp = X(:,currency:5:(end-5+currency)); % Picks only same currency
        end
    
        %% Estimate and forecast
        % Take the rolling window from the full sample
        ymean = nan(T,1);
        yfcst_lasso = nan(T,1);
        yfcst_adalasso = nan(T,1);
        yfcst_rf = nan(T,1);
        yfcst_pc = nan(T,Nfacmax);
        yfcst_spc = nan(T,Nfacmax);
        CMdltest = struct();
        for t = 1:T-troll 
            if cases == 4 
                xst_t = Xtmp(t:t+troll,:);
            else
                xst_t = X(t:t+troll,:);
            end
            yfor_t = yfor(t:t+troll);  
            
            xst_t_st = normalize(xst_t);
            % Estimate a historical mean
            ymean(troll+t) = nanmean([nan;yfor_t(1:end-2)]);
            if est_pca
           
                
                %Obtain the factor estimates
                %[Fhat,lam]=pc_factor(xst_t_st,Nf); % Yana written function
                [~,Fhat] = pca(xst_t_st,'NumComponents',Nf,'Centered',false);
                % Estimate OLS with PC factors, with varying number of components          
                y = yfor_t(1:troll);  
                for jj = 1:Nf
                  Fhat_est = Fhat(:,1:jj);
                  x = [ones(troll,1) Fhat_est(1:troll,:)];
                  beta_ols =  x'*y/size(x,1); % Since x is standardized
                  yfcst_pc(troll+t,jj) = [1 Fhat_est(end,:)]*beta_ols;  % Forecast for PC
                end
            
            end
            if est_lasso
                y = yfor_t(1:troll);
                if mod(t,re_estimate_lasso)==0 || t==1
                    try
                        if rep==1
                            rng(0, 'twister'); % This sets the seed to 0 using the Mersenne Twister generator.
                        end
                        c = cvpartition(length(y), 'KFold', 5);
                        [B,FitInfo] = lasso(xst_t_st(1:end-1,:), y, 'CV', c);
                        [x,I] = min(FitInfo.SE);
                        coeffs=B(:,I);
                        
                        
                    catch
                    end
                end
                yfcst_lasso(troll+t,1) = glmval([FitInfo.Intercept(I);coeffs],xst_t_st(end,:),'identity');
                save_coeffs_lasso(:,t).(curr) = coeffs;
            end
            if est_ada_lasso  % Adaptive LASSO
                 % Transform the predictors with respect to the weights. Use OLS
                 % estimate to calculate the weights.
                 XX = xst_t(1:end-1,:);
                 y = yfor_t(1:troll);
                 beta_OLS=(XX'*XX)\XX'*y;
                 gamma = 1;
                 w=abs(beta_OLS).^gamma;
                 w=1./w;
                 X_star=XX./w'; 
                 % Perform LASSO on the transformed data
                 if mod(t,re_estimate_lasso)==0 || t==1
                     try
                         if rep==1
                            rng(0, 'twister'); % This sets the seed to 0 using the Mersenne Twister generator.
                         end
                         c = cvpartition(length(y), 'KFold', 5);
                         [B_ada,FitInfo_ada] = lasso(X_star,y,'CV',c);
                         [x,I_ada] = min(FitInfo_ada.SE);
                         coeffs_ada=B_ada(:,I_ada)./w; %Calculate the final betas using the estimated coefficeints and the weights
                     catch
                     end
                 end
                 yfcst_adalasso(troll+t,1) = FitInfo_ada.Intercept(I_ada)+X_star(end,:)*coeffs_ada; 
                 save_coeffs_adalasso(:,t).(curr) = coeffs_ada;
   
                
            end
            if est_rf
                y = yfor_t(1:troll);
                if mod(t,re_estimate_rf)==0 || t==1
                    try
                        if rep==1
                            rng(0, 'twister'); % This sets the seed to 0 using the Mersenne Twister generator.
                        end
                        Mdl = TreeBagger(500,xst_t(1:end-1,:),y,OOBPrediction="Off",OOBPredictorImportance="Off",method="regression");
                        
                        fname = "f"+int2str(t);
                        CMdltest.(fname) = compact(Mdl); % Saves memory
                     catch
                    end
                end
            end
            if mod(t,1000) == 0
                disp("Working on observation " +t +" of " +(T-troll))
            end
           
        end % loop over t ends here
        
        if est_rf % make all forecasts for rf here (faster than in loop)
            allnames = fieldnames(CMdltest);
            tidx = [];
            for j=1:size(allnames,1)
                tidx(j) = str2double(strtok(allnames{j},'f'));
                sidx = troll+tidx(j);
                if j < size(allnames,1)
                    tidx(j+1) = str2double(strtok(allnames{j+1},'f')); 
                    eidx = troll+tidx(j+1);
                    if cases == 4 
                        yfcst_rf(sidx:eidx,1) = predict(CMdltest.(allnames{j}),Xtmp(sidx:eidx,:)); 
                    else
                        yfcst_rf(sidx:eidx,1) = predict(CMdltest.(allnames{j}),X(sidx:eidx,:)); 
                    end
                else
                    if cases == 4 
                        yfcst_rf(eidx+1:size(X,1),1) = predict(CMdltest.(allnames{j}),Xtmp(eidx+1:end,:)); 
                    else
                        yfcst_rf(eidx+1:size(X,1),1) = predict(CMdltest.(allnames{j}),X(eidx+1:end,:)); 
                    end
                end
            end
        end
  
        if est_spca %SPCAdr loops over t in the function
         
             [var_set,rmse_SPCA, Fail,Theta,r2_OS_SPCA,NR,yfcst_spc] = SPCAdr(troll,Nfacmax,yfor,X,[ones(T,1),Y(:,1)],N);           
            
             [rows, cols] = size(VAR_OUT);
             new_cols = numel([currency,var_set]);

             if new_cols > cols
                VAR_OUT = [VAR_OUT, NaN(rows, new_cols - cols)];
             end
            
             VAR_OUT(currency, 1:new_cols) = [currency, var_set];

        end
        if est_spca
            RMSE_SPCA(currency,1:size([currency,size(var_set,2),Nf,NR,Theta,rmse_SPCA*10^3],2)) = [currency,size(var_set,2),Nf,NR,Theta,rmse_SPCA*10^3]; %RMSE*10^3 measure of forecast combinations 
            R2OS_SPCA(currency,1:size([currency,size(var_set,2),Nf,NR,r2_OS_SPCA],2)) = [currency,size(var_set,2),Nf,NR,r2_OS_SPCA]; 
        end
       
        %% Calc loss functions (loss for pca is calculated in SPCAdr)
        yfcst_all = [yfcst_pc, yfcst_lasso,yfcst_adalasso, yfcst_rf, ymean];
                
        % Construct residuals %
        er = yfor - yfcst_all;
        e =er(:,1:end-1);
        e_mean =er(:,end);
        rmse_e = sqrt(nanmean(e.^2,1));
        rmse_mean = sqrt(nanmean(e_mean.^2,1));
        
        r2_OOS=100*(1-(rmse_e.^2/rmse_mean^2)); %Campbell-Thompson out-of-sample R^2
        RMSE(currency,:) = [currency,rmse_mean*10^3,rmse_e*10^3]; %RMSE*10^3 measure of forecast combinations
        R2OS(currency,:) = [currency,r2_OOS];
        
        all_forecasts.(curr) = [yfcst_all,yfcst_spc,yfor];
    end % Currency loop ends here
    %% Save results
    formatSpec = 'Case_%d';
    str = convertCharsToStrings(sprintf(formatSpec,cases));
    today = convertCharsToStrings(date);
    filename = (str)+"_"+today+"troll_"+troll+"freq_" + frequency+" samp_size="+when_end;
    %% Indicate in filenname which models are estimated
    if est_lasso
        filename = filename + "_lasso";
    end
    if est_pca
        filename = filename + "_pca";
    end
    if est_spca  
        filename = filename + "_spca";
    end
    if est_ada_lasso
        filename =filename +"_ada_lasso";
    end
    if est_rf
        filename = filename +"rf";
    end

    cd(savepath)
    save(filename,'RMSE','RMSE_SPCA','R2OS','R2OS_SPCA','VAR_OUT','all_forecasts','save_coeffs_lasso','save_coeffs_adalasso');
        
end % End of function







