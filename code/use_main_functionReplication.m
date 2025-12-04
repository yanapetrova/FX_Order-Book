% Script for 
clear; clc; close all
warning off 
[fileList, productList] = matlab.codetools.requiredFilesAndProducts('use_main_functionReplication.m');
tic

%% Pick options
reEstimateResults = 1; % If 1 re-estimates the intermediary results 
                       % This takes about 55 min for 1h frequency and 
                       % about 85 hours for the 1 minute frequency
                       % on a modern (2024) desktop with 12 physical cores

% Load data (change the path to your data location)
datalocation = 'g:\Dropbox\Yana\Data';
cd(datalocation);

% Save location  (change the path to your save data location)
savepath = 'g:\Dropbox\Yana\results';
% If you don't re-estimate the results, the supplied intermediary results
% must be placed in the savepath location.
% Note if you have several files for the same 'case' and frequency in this
% folder, the most recent will be used.

if reEstimateResults
    % This loads all data (cases are subsets of this data)
    load('datasetFX.mat')
    
    % Create 1 hour and 1 minute results (the results in tables 5 and 6)
    for freq = 1:2 % One hour and 1 minute
        if freq == 1
            frequency = 'hourly'; 
            troll = 1440; % size of rolling sample
        end
        if freq == 2
            frequency = 'minutely'; 
            troll = 7200; % size of rolling sample
        end
            for cases = 1:5 % If you use a parfor loop results with CV will not be exactly reproducable.
                if cases == 1
                    xdata = [md_IMBbest, md_IMB2, md_IMB3, md_IMB4, md_IMB5, md_dlogQask, md_dlogQbid, md_SlopeAsk, md_SlopeBid, md_Spread, md_AWS, md_OrdFlow];
                end
                if cases == 2
                    xdata = [md_IMBbest, md_IMB2, md_IMB3, md_IMB4, md_IMB5, md_dlogQask, md_dlogQbid,  md_SlopeAsk, md_SlopeBid, md_Spread, md_AWS];
                end
                if cases == 3
                    xdata =[md_IMBbest, md_IMB2, md_IMB3, md_IMB4, md_IMB5, md_dlogQask, md_dlogQbid, md_SlopeAsk, md_SlopeBid, md_Spread];
                end
                if cases == 4 % will pick data for currency later in the function
                    xdata =[md_IMBbest, md_IMB2, md_IMB3, md_IMB4, md_IMB5, md_dlogQask, md_dlogQbid, md_SlopeAsk, md_SlopeBid, md_Spread, md_AWS, md_OrdFlow];
                end
                if cases == 5 % will pick data for currency later in the function for case 5
                    xdata = md_OrdFlow;
                end
                
                disp("Working on case: "+cases +" freq: " +frequency +" Troll: "+troll)
                main_function_replication(xdata,md_Y,cases,troll,frequency,savepath)
            end
    end
 
    toc
end
% Build the tables from the intermediate results
% The location of the files are set in savepath above
% Currencies
curr = {'USD/CAD'; 'AUD/USD'; 'GBP/USD'; 'USD/JPY'; 'EUR/USD'};

% Column names for numeric columns
cols = {'pc1','pc2','pc3','pc4','pc5', ...
        's_pc1','s_pc2','s_pc3','s_pc4','s_pc5', ...
        'LASSO','ada_LASSO','RF'};

nCases = 5;
nCurr  = numel(curr);
nCols  = numel(cols);

% Preallocate row container
rows_Label = strings(nCases*(nCurr+1),1);   % +1 header row per case
rows_Data  = nan(nCases*(nCurr+1), nCols);

idx = 1;

for c = 1:nCases

    % ---- Case header row ----
    rows_Label(idx) = sprintf("Case %d", c);
    rows_Data(idx,:) = nan(1,nCols);
    idx = idx + 1;

    % ---- Currency rows ----
    for k = 1:nCurr
        rows_Label(idx) = curr{k};
        rows_Data(idx,:) = nan(1,nCols);
        idx = idx + 1;
    end
end

% Final table
ResultTable = array2table(rows_Data, 'VariableNames', cols);
ResultTable.Label = rows_Label;

% Move Label to front
ResultTable5 = movevars(ResultTable, 'Label', 'Before', 1);
ResultTable6 = ResultTable5;
cases = [1:5];

latestMinutely = cell(numel(cases),1);
latestHourly   = cell(numel(cases),1);
for k = 1:numel(cases)
    c = cases(k);

    patHour = sprintf('*case_%d*hourly*', c);
    Lhour   = dir(fullfile(savepath, patHour));
    if ~isempty(Lhour)
        % pick the most recently modified
        [~, idxHour] = max([Lhour.datenum]);
        load(fullfile(savepath, Lhour(idxHour).name));
        ResultTable5{k*5-4+k:k*5+k,2:6} = R2OS(:,2:6);
        ResultTable5{k*5-4+k:k*5+k,7:11} = R2OS_SPCA(:,6:10);
        ResultTable5{k*5-4+k:k*5+k,12:14} = R2OS(:,7:9);
        
    else
        latestHourly{k} = '';
        warning('No hourly files found for case_%d.', c);
    end

    patMin = sprintf('*case_%d*minutely*', c);
    Lmin   = dir(fullfile(savepath, patMin));     

    if ~isempty(Lmin)
        % pick the most recently modified
        [~, idxMin] = max([Lmin.datenum]);
        latestMinutely{k} = fullfile(savepath, Lmin(idxMin).name);
        load(fullfile(savepath, Lmin(idxMin).name));
        ResultTable6{k*5-4+k:k*5+k,2:6} = R2OS(:,2:6);
        ResultTable6{k*5-4+k:k*5+k,7:11} = R2OS_SPCA(:,6:10);
        ResultTable6{k*5-4+k:k*5+k,12:14} = R2OS(:,7:9);
    else
        latestMinutely{k} = '';
        warning('No minutely files found for case_%d.', c);
    end
end

%% Construct figures 3 and 4, and table 7
c = 1; % First case
patMin = sprintf('*case_%d*minutely*', c);
Lmin   = dir(fullfile(savepath, patMin));
[~, idxMin] = max([Lmin.datenum]);

% Lasso 1 minute
results = load(fullfile(savepath, Lmin(idxMin).name));
varnames = [repmat({'IMBbest'},1,5), repmat({'IMB2'},1,5)  repmat({'IMB3'},1,5), repmat({'IMB4'},1,5), ...
    repmat({'IMB5'},1,5), repmat({'dlogQask'},1,5), repmat({'dlogQbid'},1,5), repmat({'SlopeAsk'},1,5), repmat({'SlopeBid'},1,5), ...
    repmat({'Spread'},1,5), repmat({'AWS'},1,5), repmat({'OrdFlow'},1,5) ];
currencies = [{'USDCAD'},{'AUDUSD'},{'GBPUSD'},{'USDJPY'},{'EURUSD'}];
for k=1:12
    varnames(2,1+5*(k-1):5+5*(k-1)) = currencies;
end

% Initialize an empty 1x60 cell array
names = cell(1, 60);

% Loop through each column and concatenate the strings vertically
for i = 1:60
    tmp = strjoin(varnames(:, i));
    names{i} = strrep(tmp,' ','-');
end
varnames = [repmat({'IMBbest'},1,5), repmat({'IMB2'},1,5)  repmat({'IMB3'},1,5), repmat({'IMB4'},1,5), ...
    repmat({'IMB5'},1,5), repmat({'dlogQask'},1,5), repmat({'dlogQbid'},1,5), repmat({'SlopeAsk'},1,5), repmat({'SlopeBid'},1,5), ...
    repmat({'Spread'},1,5), repmat({'AWS'},1,5), repmat({'OrdFlow'},1,5) ];
currencies = [{'USDCAD'},{'AUDUSD'},{'GBPUSD'},{'USDJPY'},{'EURUSD'}];
for k=1:12
    varnames(2,1+5*(k-1):5+5*(k-1)) = currencies;
end

% Initialize an empty 1x60 cell array
names = cell(1, 60);

% Loop through each column and concatenate the strings vertically
for i = 1:60
    tmp = strjoin(varnames(:, i));
    names{i} = strrep(tmp,' ','-');
end

% Construct out of sample R2 for each nonoverlapping window EURUSD
Actual = results.all_forecasts.EURUSD(:,end);
Lasso = results.all_forecasts.EURUSD(:,6); 
meanFor = results.all_forecasts.EURUSD(:,9);
i = 1;
for j=7201:7200:(size(Actual,1)-7200)
    e = Actual(j:(j+7199)) - Lasso(j:(j+7199));
    e_mean = Actual(j:(j+7199)) - meanFor(j:(j+7199));
    rmse_e = sqrt(nanmean(e.^2,1));
    rmse_mean = sqrt(nanmean(e_mean.^2,1));
    r2_OOS(i,:) = 100*(1-(rmse_e.^2/rmse_mean^2)); %Campbell-Thompson out-of-sample R^2
    i =  i + 1;
end
r2_OOSEURUSD = r2_OOS;


% Construct out of sample R2 for each nonoverlapping window AUDUSD
Actual = results.all_forecasts.AUDUSD(:,end);
Lasso = results.all_forecasts.AUDUSD(:,6); 
meanFor = results.all_forecasts.AUDUSD(:,9);
i = 1;
for j=7201:7200:(size(Actual,1)-7200)
    e = Actual(j:(j+7199)) - Lasso(j:(j+7199));
    e_mean = Actual(j:(j+7199)) - meanFor(j:(j+7199));
    rmse_e = sqrt(nanmean(e.^2,1));
    rmse_mean = sqrt(nanmean(e_mean.^2,1));
    r2_OOS(i,:) = 100*(1-(rmse_e.^2/rmse_mean^2)); %Campbell-Thompson out-of-sample R^2
    i =  i + 1;
end
r2_OOSAUDUSD = r2_OOS;

save_coeffs_lasso = results.save_coeffs_lasso;
save_coeffs_adalasso = results.save_coeffs_adalasso;
fieldNames = fieldnames(save_coeffs_lasso);
for i = 1:5 % 
    currentField = fieldNames{i};
    % Graphs
    for j=1:size(save_coeffs_lasso,2)
        fx(:,j) = save_coeffs_lasso(j).(currentField);
        fxada(:,j) = save_coeffs_adalasso(j).(currentField);
    end
    fx(fx~=0) = 1;
    fxada(fxada~=0) = 1;
    %% Pick every 7200th column to get non-overlapping windows
    p = 0;
    for k=1:7200:size(fx,2) 
        p = p + 1;
        nonoverlappingUSDCAD(:,p) = fx(:,k);
    end
     
    % Convert start and end dates to datetime format
    startDate = datetime('2019-10-08', 'Format', 'yyyy-MM-dd');
    endDate = datetime('2020-04-19', 'Format', 'yyyy-MM-dd');
    
    % Calculate the time interval between start and end dates
    totalDays = days(endDate - startDate);
    % Divide the interval into 7 equal parts to get 8 equidistant points
    intervalDays = totalDays / 7;
    % Generate the 8 equidistant dates
    equidistantDates = startDate + days(0:intervalDays:totalDays);
    % Convert to desired format
    equidistantDatesStr = datestr(equidistantDates, 'yyyy/mm/dd');
    % find the most kept coeffecs for plotting
    [kept_often,I] = maxk(sum(fx,2),10);
    [kept_oftenada,Iada] = maxk(sum(fxada,2),10);
    lassotable1min((3*i)-2:3*i,1:3) = [repmat({currentField},3,1),names(1,I(1:3))',num2cell(100*kept_often(1:3)./size(fx,2))];
    lassotable1min((3*i)-2:3*i,4:5) = [names(1,Iada(1:3))',num2cell(100*kept_oftenada(1:3)./size(fx,2))];
    
    
    if i == 2 || i == 5 % Only plot two figures for paper
        fig = figure;
        set(fig, 'Position', [100, 100, 900, 450]);
        spy(nonoverlappingUSDCAD(I,:));
        hold on
        
        nRows = numel(I);           % e.g. 10
        nW    = size(nonoverlappingUSDCAD,2);  % e.g. 26
        
        % 1) place R^2 numbers ABOVE the first row:
        yText = 0.5;  % half‐unit *above* row 1
        if i == 2
            r2_OOS = r2_OOSAUDUSD;
        end
        if i == 5
            r2_OOS = r2_OOSEURUSD;
        end
    
        r2labels = arrayfun(@(x) sprintf('%.1f',x), r2_OOS, 'UniformOutput',false);
        text(1:nW-1, repmat(yText,1,nW-1), r2labels, ...
             'HorizontalAlignment','center', ...
             'VerticalAlignment','middle', ...
             'FontSize',9);
        
        % 2) stick the “R²” label just to the left of column 1
        text(0.5, yText, '$R^2$', ...
             'HorizontalAlignment','right', ...
             'VerticalAlignment','middle', ...
             'Interpreter','latex', ...
             'FontSize',10);
        
        % 3) expand the axes so 0.5 isn’t clipped
        xlim([0.5, nW+0.5]);
        ylim([0,   nRows+0.5]);    % note lower bound is 0 now
        
        hold off
        
        % 4) now re-set your ticks & labels
        yticks(1:nRows);
        yticklabels(names(I));              
          
        xticks(linspace(1,nW,8));
        xticklabels(equidistantDatesStr);  
    
        if i == 5
            title("EURUSD (fig 4)")
            saveas(fig, savepath+"\fig4.fig");
        end
        if i == 2
            title("AUDUSD (fig 3)")
            saveas(fig, savepath+"\fig3.fig");
        end
    end

end
ResultTable7 = cell2table(lassotable1min);
ResultTable7.Properties.VariableNames = {'Currencies','Lasso','%incl','ada-LASSO','%includ'};

writetable(ResultTable5,savepath+"\Table5.xlsx")
writetable(ResultTable6,savepath+"\Table6.xlsx")
writetable(ResultTable7,savepath+"\Table7.xlsx")



