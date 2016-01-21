classdef mPortfolio < handle
    % mPORTFOLIO Portfolio(投资组合)类
    %   Detailed explanation goes here
    
    % define properties
    % 定义属性
    % equity: 权益的时间序列
    % positions: 仓位数据表，包含所有交易对象的仓位信息
    % trades: 逐笔交易数据表，包含所有交易对象的交易信息
    properties
        equity
        returns
        PosData % 所有交易的金融工具的仓位信息数据
        TradeData % 所有交易的金融工具的交易成交数据
        Instruments %
        positions % 仓位数据，每个金融工具分配单独的一个table
        trades % 成交数据，每个金融工具分配单独的一个table
        drawdownStats % 回撤数据
        perTradeStats % 储存从开仓到平仓的逐笔交易数据统计
        txnPL % 储存每笔成交数据的盈亏计算
        positionsPL %储存每个 symbol 的仓位盈亏计算结果
        symbolPL % 储存每个 Symbol 上的合计盈亏
        portfolioPL % 储存投资组合层面的盈亏数据
        performance % 储存业绩分析的各种 ratio
        basicStats % 收益率的基本统计量
        
    end
    
    methods
        % constructor
        function portf = mPortfolio(equity,posData,tradeData)
            % 检查输入的posData和tradeData中的symbol是否一致
            % 记录仓位和成交记录的数据一致性的第一步是检查所包含的金融工具的名字（即Symbol字段）是否完全一致
            % 仓位数据是可以通过成交明细计算出的，所以从成交明细计算得出的仓位数据可以和posData比较以检查posData和tradeData的一致性
            % 此处只进行Symbol字段的检查
            % 如果tradeData和posData中的symbol个数不一致，说明数据可能有错漏，或者tradeData中包含了日内交易数据
            if length( unique( posData.Symbol ) ) ~= length( unique( tradeData.Symbol ) )
            warning('posData和tradeData中的symbol个数不一致,可能tradeData存在日内交易数据,而posData为每日结算后的仓位')
            end
            
            % 属性的初始化
            portf.equity = equity; %每日权益
            portf.PosData = posData;
            portf.TradeData = tradeData;
            
            portf.Instruments = struct(); % 用于储存每个instrument的属性如 Multiplier
            portf.positions = struct(); % 用于储存每个instrument的position data
            portf.trades = struct(); % 用于储存每个instrument的trade data
            portf.performance = struct(); % 用于储存业绩分析的结果
            portf.basicStats = struct(); % 用于储存收益序列的基础统计量
            portf.returns = struct();
            portf.drawdownStats = table();
            portf.perTradeStats = struct();
            portf.txnPL = struct();
            portf.symbolPL = struct();
            portf.portfolioPL = struct();
            portf.positionsPL = struct();
            
            disp( '成功创建mPortfolio类对象' )
        end
        
        % 向portfolio添加所交易的金融工具的信息数据，从tradeData中生成
        function addInstruments(portfolio, varargin)
            if isempty(portfolio.TradeData)
                error('TradeData is empty, cannot create instruments from TradeData. TradeData是空的')
            end
            if nargin == 1
                error('Please input the variable we need in TradeData. 请输入字段名') 
            end
            
            labels = { 'Symbol',char(varargin)};
            [Lia,~] = ismember( labels,fieldnames(portfolio.TradeData) );
            if ~all(Lia)   
                error('输入的值中有的不是TradeData的字段')
            end
            temp = unique(portfolio.TradeData(:,labels) );
            symbolnames = unique(temp.Symbol);
            if size(temp,1) ~= length( symbolnames )
                error('存在金融工具有不同的属性数据，请检查TradeData')
            end
            for i=1:length(symbolnames)
               
               portfolio.addInstrument(temp.Symbol(i),temp(i,varargin) );
            end
            
            disp( ['生成了',num2str(length(symbolnames)),'个Instrument数据'] )
        end % end of function addInstruments
        
        
        % 向portfolio添加所交易的金融工具的信息数据，根据自定义的数据添加
        function addInstrument(portfolio, symbol, da)
            symbol = char(symbol);
            portfolio.Instruments.(symbol) = da;
        end
        
        % create positions
        % 创建 positions 中的仓位数据表，每个instrument分配一个单独的表
        % 从 PosData 中按symbol分割数据表
        function cpositions(portfolio)
            symbolnames = unique( portfolio.PosData.Symbol );
            for i = 1:length(symbolnames)
               symbol = symbolnames(i);
              %idx = find(portfolio.PosData.Symbol==symbol);
               portfolio.positions.(char(symbol)) = portfolio.PosData(portfolio.PosData.Symbol==symbol,:);
            end
            
            disp( ['生成了',num2str(length(symbolnames)),'个position数据表'] )
        end % end of cpositions
        
        % create trades
        % 创建 trades 中的成交记录数据表，每个instrument分配一个单独的表
        % 从 TradeData 中按symbol分割数据表
        function ctrades(portfolio)
            symbolnames = unique( portfolio.TradeData.Symbol );
            for i = 1:length(symbolnames)
               symbol = symbolnames(i);
              %idx = find(portfolio.PosData.Symbol==symbol);
               portfolio.trades.(char(symbol)) = portfolio.TradeData(portfolio.TradeData.Symbol==symbol,:);
            end
            
            disp( ['生成了',num2str(length(symbolnames)),'个trade数据表'] )
        end % end of ctrades
        
         % 基本统计量
        function calcBasicStats(portfolio,period)
            % period的默认值是daily
            if nargin == 1
               period = 'daily'; 
            end
            if ~ismember(period,{'daily','weekly','monthly'})
                error('period的输入值只能是 daily,weekly,monthly')
            end
            period=char(period);
            % mean
            rtnMean = ...
                mean(fts2mat(portfolio.returns.(period)));
            
            % median
            rtnMedian = ...
                median( fts2mat(portfolio.returns.(period)));
            
            % std
            rtnStd = ...
                std( fts2mat(portfolio.returns.(period)));
            annRtn = portfolio.getAnnReturn(period);
            portfolio.basicStats.(period) = table( rtnMean,rtnMedian,rtnStd,annRtn );
            disp( [period,' return的basicStats：'])
            disp(portfolio.basicStats.(period))
            
        end % end of calcBasicStats
        
        % 计算简单收益率；根据equity变量中的equity序列的权益值计算，周期等于权益值序列的周期
        % returns序列的第一个值为 NaN
        function calcDailyReturn(portfolio,method)
            % method的默认值是Simple
            if nargin == 1
               method = 'Simple'; 
            end
             % 检查日权益序列是否为空
            if isempty(portfolio.equity)
               error('equity is empty. equity为空.') 
            end
            % 检查 method 的值是否为Simple或Continuous
            if ~ismember(method,{'Simple','Continuous'})
                error('method的输入值只能是 Simple 或 Continuous')
            end
             portfolio.returns.daily = tick2ret(portfolio.equity.equity,'method',method) ;
             
        end
        
        function calcWeeklyReturn(portfolio,method)
            % method的默认值是Simple
            if nargin == 1
               method = 'Simple'; 
            end
            % 检查日权益序列是否为空
            if isempty(portfolio.equity)
               error('equity is empty. equity为空.') 
            end
            % 检查 method 的值是否为Simple或Continuous
            if ~ismember(method,{'Simple','Continuous'})
                error('method的输入值只能是 Simple 或 Continuous')
            end
            % 将日权益序列转化为周度数据
            equity_weekly = toweekly(portfolio.equity);
            % 计算月收益率
            portfolio.returns.weekly = tick2ret(equity_weekly,'method',method);
        end
        
        function calcMonthlyReturn(portfolio,method)
            % method的默认值是Simple
            if nargin == 1
               method = 'Simple'; 
            end
            % 检查日权益序列是否为空
            if isempty(portfolio.equity)
               error('equity is empty. equity为空.') 
            end
            % 检查 method 的值是否为Simple或Continuous
            if ~ismember(method,{'Simple','Continuous'})
                error('method的输入值只能是 Simple 或 Continuous')
            end
            % 将日权益序列转化为月度数据
            equity_monthly = tomonthly(portfolio.equity);
            % 计算月收益率
            portfolio.returns.monthly = tick2ret(equity_monthly,'method',method);
        end % end of calcMonthlyReturn
        
        % 计算夏普比率 sharpe ratio
        % riskfree 为无风险利率, default 0
        % period 为收益率的周期, default 'daily'
        function [shpr] = calcSharpe(portfolio,riskfree,period)
            if nargin == 1
               riskfree = 0; % risk free rate default to be zero
               period = 'daily'; % period default to be daily 
            end
            if nargin == 2
               period = 'daily'; % period default to be daily 
            end
            if ~ismember(period,{'daily','weekly','monthly','yearly'})
                error('period的输入值只能是 daily,weekly,monthly,yearly')
            end
            
            % 对 risk free rate 转换，使周期和 return 的周期一致
            T = 250;
            switch period
                case 'daily'
                    T = 250;
                    
                case 'weekly'
                    T = 52;
                    
                case 'monthly'
                    T = 12;
                case 'yearly'
                    T =1;
            end
            riskfree = (1+riskfree)^(1/T) - 1 ; 
            portfolio.performance.(char(period)).sharpe = sharpe(fts2mat(portfolio.returns.(char(period))), riskfree);
            
            % 根据输出变量选择输出方式，如引用函数时语句有输出变量，给输出变量赋值，否则在控制台打印计算结果
            if nargout == 1
            shpr = portfolio.performance.(char(period)).sharpe;
            else 
            disp( [ period,' return的sharpe ratio为', num2str(portfolio.performance.(char(period)).sharpe) ] )
            end
        end
        
        function calcMaxDrawdown(portfolio)
            portfolio.performance.daily.maxDrawdown = maxdrawdown(fts2mat(portfolio.equity.equity),'return');
            disp( [ 'MaxDrawdown为', num2str(portfolio.performance.daily.maxDrawdown) ] )
        end
        
        % 计算 sortino ratio
        % MAR: 年化的目标收益率，投资者常将其设为 0, 历史均值收益率, 无风险利率, 或基准利率
        % 调用下跌半偏差函数 downsideDeviation(x,target)
        % target 为将MAR调整为与x的周期对应的值. 比如：x 是月收益率，则target为调整后的月目标收益率
        function [sortino] = calcSortino(portfolio, MAR, period)
            % set default values, 设置默认值
            if nargin == 1
               MAR = 0;
               period = 'monthly';
            end
            if nargin == 2
               period = 'monthly'; 
            end
            % 调整目标收益率
            T = 12;
            switch period
                case 'daily'
                    T = 250;
                    
                case 'weekly'
                    T = 52;
                    
                case 'monthly'
                    T = 12;
                    
            end
            target = ( 1+MAR )^( 1/T ) - 1 ;
            annRtn = portfolio.getAnnReturn(period);
            dd = downsideDeviation( fts2mat(portfolio.returns.(char(period))), target );
            annDd = sqrt(T)*dd;
            portfolio.performance.(char(period)).sortino = (annRtn - MAR)/annDd;
            sortino = (annRtn - MAR)/annDd;
            if nargout == 0
                disp( [ period,' return的以',num2str(MAR),'为目标收益率的Sortino ratio为', num2str(portfolio.performance.(char(period)).sortino) ] )
            end
        end % end of calcSortino
        
        % upside potential ratio
        % 有些投资者寻找的不是一段时期内平均收益最高的经理，而是那些高于MAR的平均收益率最高的经理
        % Sortino, van der Meer and Plantinga(1999)
        % 提出了用上行空间来代替Sortino比率中的分子中的额外收益
        % 上行空间指的是超过MAR的预期收益率，可以看做是成功的潜力
        % upside potential ratio = E[ r - MAR ]/downside risk
        % 和 Sortino 一样，upside potential ratio 也是年化的数据
        function [uppr] = calcUpsidePotential(portfolio,MAR,period)
             % set default values, 设置默认值
            if nargin == 1
               MAR = 0;
               period = 'monthly';
            end
            if nargin == 2
               period = 'monthly'; 
            end
             % 调整目标收益率
            T = 12;
            switch period
                case 'daily'
                    T = 250;
                    
                case 'weekly'
                    T = 52;
                    
                case 'monthly'
                    T = 12;
                    
            end
            target = ( 1+MAR )^( 1/T ) - 1 ;
            dd = downsideDeviation( fts2mat(portfolio.returns.(char(period))), target );
            excess = fts2mat(portfolio.returns.(char(period))) - target;
            uppr = mean(max(excess,0)) / dd;
            portfolio.performance.(char(period)).upsidePotential = uppr ;
            
            if nargout == 0
                disp( [ 'Upside potential ratio 基于 ',period,' return： ',num2str(uppr) ] )
            end
        end
        
        
        % Sterling ratio 计算Sterling比率
        % Sterling和Burke比率受到CTAs的广泛欢迎，因为这两个比率呈现了他们所相信做的最好的那方面：
        % 即利润最大化，但必须严格控制损失.
        % Sterling比率考虑用回撤来度量风险，所以比Sortino比率更进了一步
        % sterling = (mean(r) - riskfree )/mean(drawdown)
        % 分母被定义为数值显著的一些回撤的平均值，"数值显著"还有待定义.此处将其定义为所有回撤的平均值
        function [sterling] = calcSterling(portfolio,riskfree,period)
            
            if nargin == 1
               riskfree = 0; % risk free rate default to be zero
               period = 'monthly'; % period default to be daily 
            end
            if nargin == 2
               period = 'monthly'; % period default to be daily 
            end
            if ~ismember(period,{'daily','weekly','monthly','yearly'})
                error('period的输入值只能是 daily,weekly,monthly,yearly')
            end
         
            eqty = fts2mat(portfolio.equity);
            eqty = eqty(:,1);
            % cumulative maximum
            % index of cmmax is the same as eqty
            cmmax = cummax(eqty,1); 
            % 回撤计算公式： drawdown :=  cummax - eqty (which is positive)
            % dd and ddrate have same indexing as eqty
            dd =  cmmax - eqty ;
            ddrate = dd ./cmmax ;
            
            
            annRtn = portfolio.getAnnReturn( period );
            mdwn = mean( ddrate );
            sterling = ( annRtn - riskfree)/mdwn;
            portfolio.performance.(char(period)).sterling = sterling;
        end % end of calcSterling
        
        
        %用各个回撤期的最大回撤代替所有的回撤
        function [sterling_maxDD] = calcSterling_maxDD(portfolio,riskfree,period)
            
            if nargin == 1
               riskfree = 0; % risk free rate default to be zero
               period = 'monthly'; % period default to be daily 
            end
            if nargin == 2
               period = 'monthly'; % period default to be daily 
            end
            if ~ismember(period,{'daily','weekly','monthly','yearly'})
                error('period的输入值只能是 daily,weekly,monthly,yearly')
            end
  
            annRtn = portfolio.getAnnReturn(period);
            mdwn = mean( portfolio.drawdownStats.maxDrawdownPct./100 );
            sterling_maxDD = ( annRtn - riskfree)/mdwn;
            portfolio.performance.(char(period)).sterling_maxDD = sterling_maxDD;
        end % end of calcSterling_maxDD
        
        % Burke ratio
        % Burke类似于Sterling ratio,
        % Burke(1994)建议使用每个回撤的平方和的平方根以惩罚与大量温和回撤相对的幅度较大的回撤
        function [ burke ] = calcBurke(portfolio, riskfree,period)
            if nargin == 1
               riskfree = 0; % risk free rate default to be zero
               period = 'monthly'; % period default to be daily 
            end
            if nargin == 2
               period = 'monthly'; % period default to be daily 
            end
            if ~ismember(period,{'daily','weekly','monthly','yearly'})
                error('period的输入值只能是 daily,weekly,monthly,yearly')
            end
          
            eqty = fts2mat(portfolio.equity);
            eqty = eqty(:,1);
            % cumulative maximum
            % index of cmmax is the same as eqty
            cmmax = cummax(eqty,1); 
            % 回撤计算公式： drawdown :=  cummax - eqty
            % dd and ddrate have same indexing as eqty
            dd =  cmmax - eqty ;
            ddrate = dd ./cmmax ;
            
          
            deno = sqrt( sum(ddrate.^2) );
            burke = (portfolio.getAnnReturn(period) - riskfree)/deno;
            portfolio.performance.(char(period)).burke = burke;
        end % end of calcBurke        
        
        
        
        function [ burke_maxDD ] = calcBurke_maxDD(portfolio, riskfree, period)
            if nargin == 1
               riskfree = 0; % risk free rate default to be zero
               period = 'monthly'; % period default to be daily 
            end
            if nargin == 2
               period = 'monthly'; % period default to be daily 
            end
            if ~ismember(period,{'daily','weekly','monthly','yearly'})
                error('period的输入值只能是 daily,weekly,monthly,yearly')
            end
         
            dwn = portfolio.drawdownStats.maxDrawdownPct./100 ;
            deno = sqrt( sum(dwn.^2) );
            burke_maxDD = (portfolio.getAnnReturn(period) - riskfree)/deno;
            portfolio.performance.(char(period)).burke_maxDD = burke_maxDD;
        end % end of calcBurke_maxDD
        
        
        
        % 统计每一个回撤期的数据
        function calcDrawdownStats(portfolio)
            % 检查日权益序列是否为空
            if isempty(portfolio.equity)
               error('equity is empty. equity为空.') 
            end
            % convert fts to array
            eqty = fts2mat(portfolio.equity);
            eqty = eqty(:,1);
            % cumulative maximum
            % index of cmmax is the same as eqty
            cmmax = cummax(eqty,1); 
            % 回撤计算公式： drawdown :=  cummax - eqty
            % dd and ddrate have same indexing as eqty
            dd =  cmmax - eqty ;
            ddrate = dd ./cmmax ;
            % 回撤日的index. get the index number of drawdown days
            % for indexing in eqty
            ddidx = find(dd>0);
          
            startidx = [ddidx(1);ddidx( [1;diff(ddidx)]>1 )];
            endidx = [ddidx( [diff(ddidx);1]>1 ) ; ddidx(end)] ;
            startdate = cellstr(datestr(portfolio.equity.dates(startidx),'yyyy-mm-dd'));
            enddate =  cellstr(datestr(portfolio.equity.dates(endidx),'yyyy-mm-dd'));
            timespan = endidx - startidx +1;
            %-% 每个回撤期内可能有多个局部高点，寻找local maximum的方法不适用
            % [maxDrawdown, maxDDloc] = findpeaks(-dd);
            % 
            %maxDrawdownDate = datestr(portfolio.equity.dates(maxDDloc));
            %maxDrawdownRate = ddrate(maxDDloc);
            %-%
            nrows = length(startidx);
            maxDrawdown = zeros(nrows,1); % double
            maxDrawdownPct = zeros(nrows,1); % double
            maxDrawdownDate = zeros(nrows,1); % datenum
            for i =1:nrows
               dd_temp = dd(startidx(i):endidx(i)); 
               
               [maxDrawdown(i,1),idx_temp] = max(dd_temp);
               maxDrawdownDate(i,1) = portfolio.equity.dates( startidx(i)+idx_temp-1);
               maxDrawdownPct(i,1) = 100*ddrate(startidx(i)+idx_temp-1);
            end
            maxDrawdownDate = cellstr(datestr(maxDrawdownDate,'yyyy-mm-dd')) ;
            portfolio.drawdownStats = table(startdate,enddate,...
                maxDrawdown,maxDrawdownPct,maxDrawdownDate,timespan);
            disp('计算回撤数据完成,可调用drawdownStats属性查看....')
        end
        
        % 计算下跌频率
        function [dwnfreq] = calcDownsideFreq(portfolio,MAR,period)
            if nargin ==1
               MAR = 0;
               period = 'monthly';
            end
            if nargin == 2
               period = 'monthly'; 
            end
            if ~ismember(period,{'daily','weekly','monthly','yearly'})
                error('period的输入值只能是 daily,weekly,monthly,yearly')
            end
            
             % 调整目标收益率
            T = 12;
            switch period
                case 'daily'
                    T = 250;
                    
                case 'weekly'
                    T = 52;
                    
                case 'monthly'
                    T = 12;         
            end
            target = ( 1+MAR )^( 1/T ) - 1 ;
            rtn = fts2mat( portfolio.returns.(char(period)) );
            rtn_down = rtn( rtn < target );
            dwnfreq = length( rtn_down )/length(rtn) ;
            if nargout == 0
               disp( [period,' return Downside Frequency is: ',num2str(dwnfreq)] ) 
            end
        end
        
        % 计算 正收益的sharpe
        function [sharpe_upside] = calcSharpe_upside(portfolio,MAR,period)
            
              % set default values, 设置默认值
            if nargin == 1
               MAR = 0;
               period = 'monthly';
            end
            if nargin == 2
               period = 'monthly'; 
            end
             % 调整目标收益率
            T = 12;
            switch period
                case 'daily'
                    T = 250;
                    
                case 'weekly'
                    T = 52;
                    
                case 'monthly'
                    T = 12;
                    
            end
            target = ( 1+MAR )^( 1/T ) - 1 ;
            
            rtn = fts2mat(portfolio.returns.(char(period)));
            rtn_posi = rtn( rtn>target );
            sharpe_upside = sharpe(rtn_posi);
            
            portfolio.performance.(char(period)).sharpe_upside = sharpe_upside;
        end
        
        function [sharpe_downside] = calcSharpe_downside(portfolio,MAR, period)
                          % set default values, 设置默认值
            if nargin == 1
               MAR = 0;
               period = 'monthly';
            end
            if nargin == 2
               period = 'monthly'; 
            end
             % 调整目标收益率
            T = 12;
            switch period
                case 'daily'
                    T = 250;
                    
                case 'weekly'
                    T = 52;
                    
                case 'monthly'
                    T = 12;
                    
            end
            target = ( 1+MAR )^( 1/T ) - 1 ;
            
            rtn = fts2mat(portfolio.returns.(char(period)));
            rtn_nega = rtn( rtn<=target );
            sharpe_downside = sharpe(rtn_nega);
            
            portfolio.performance.(char(period)).sharpe_downside = sharpe_downside;
        end
        
        function calcTxn_allInstruments(portfolio)
            symbols = fieldnames(portfolio.Instruments);
            nsymbols = length(symbols);
            for i=1:nsymbols
               portfolio.calcTxn_s(symbols(i)); 
            end
            disp(['完成了 ',num2str(nsymbols),'个对象的交易数据的计算'])
        end
        
        function [ result ] = calcTxn_s(portfolio,symbol)
            n = size(portfolio.trades.(char(symbol)),1);
            portfolio.trades.(char(symbol)).txnId = (1:n)' ;
            % 计算多头
            txn_long = portfolio.calcTxn_long(symbol);
            
            % 计算空头
            txn_short = portfolio.calcTxn_short(symbol);
            % 合并多空
            result = [txn_long;txn_short];
            result = sortrows(result,'txnId','ascend');
            % 计算 多空持仓
            longPos = cumsum(result.LongQty);
            longPos = table(longPos);
            shortPos = cumsum(result.ShortQty);
            shortPos = table(shortPos);
            result = [result,longPos,shortPos];
            % 计算 cumPL
            cumPL = cumsum( result.realizedPL);
            cumPL = table(cumPL);
            
            cumNetPL = cumsum( result.netPL);
            cumNetPL = table(cumNetPL);
            result = [result,cumPL,cumNetPL];
            portfolio.txnPL.(char(symbol)) = result ;
        end
        
        function [ res ] = calcTxn_long(portfolio,symbol)
            tradeda = portfolio.trades.(char(symbol)) ;
            tradeda = tradeda( tradeda.LongQty ~= 0,: ) ;
            nrows = size(tradeda, 1);
            if nrows == 0
                res = table();
                return
            end
            action = zeros(nrows,1);
            action = categorical(cellstr(char(action))) ;
            
            
            realizedPL = zeros(nrows,1);
            prevPos=zeros(nrows,1);
            
            avgPosCost = zeros(nrows,1);
            
            % 计算 current position and previous position
            curPos = cumsum(tradeda.LongQty);
            % previous position of first record is default to be 0, unless
            % we need to consider situation that the portfolio has position
            % before the performence analysis period. and this will be
            % continued to develop in the future.
            prevPos(1) = 0;
            prevPos(2:end) = curPos(1:end-1);
            
            % 计算 action
            action(prevPos==0)='open';
            action( curPos==0 ) ='flatten';
            action( prevPos~=0 & curPos ~=0 & curPos > prevPos)='increase';
            action( prevPos~=0 & curPos ~=0 & curPos < prevPos)='decrease';
            
            % 计算每笔交易之后的平均持仓成本 
            avgPosCost(1) = tradeda.LongPrice(1);
            for i=2:nrows
               switch action(i)
                   case 'open'
                       avgPosCost(i) = tradeda.LongPrice(i);
                        
                   case 'increase'
                       avgPosCost(i) = ...
                          ( avgPosCost(i-1)*prevPos(i)+ tradeda.LongPrice(i)*tradeda.LongQty(i))/curPos(i);
                   case 'decrease'
                       avgPosCost(i) = avgPosCost(i-1);
                   case 'flatten'
                       avgPosCost(i) = avgPosCost(i-1);
               end
            end
            
            % 计算 realizedPL
            temp_idx = find(action=='decrease'|action=='flatten');
            realizedPL(temp_idx) =  (tradeda.LongPrice(temp_idx) - avgPosCost(temp_idx-1)).* abs(tradeda.LongQty(temp_idx))...
                *portfolio.Instruments.(char(symbol)).Multiplier ;
            netPL = realizedPL - tradeda.TxnFee;
                
            res = [tradeda,table(action,prevPos,curPos,avgPosCost,netPL,realizedPL)];
        end
        
        function [res] = calcTxn_short(portfolio,symbol)
            tradeda = portfolio.trades.(char(symbol)) ;
            tradeda = tradeda( tradeda.ShortQty ~= 0,: ) ;
            nrows = size(tradeda, 1);
            if nrows == 0
                res = table();
                return
            end
            action = zeros(nrows,1);
            action = categorical(cellstr(char(action))) ;
            
            
            realizedPL = zeros(nrows,1);
            prevPos=zeros(nrows,1);
            
            avgPosCost = zeros(nrows,1);
            
            
            curPos = cumsum(tradeda.ShortQty);
            prevPos(1) = 0;
            prevPos(2:end) = curPos(1:end-1);
            
            % 计算 action
            action(prevPos==0) = 'open';
            action( curPos==0 ) = 'flatten';
            action( prevPos~=0 & curPos ~=0 & curPos > prevPos)='decrease';
            action( prevPos~=0 & curPos ~=0 & curPos < prevPos)='increase';
            
            % 计算每笔交易之后的平均持仓成本 
            avgPosCost(1) = tradeda.ShortPrice(1);
            for i=2:nrows
               switch action(i)
                   case 'open'
                       avgPosCost(i) = tradeda.ShortPrice(i);
                        
                   case 'increase'
                       avgPosCost(i) = ...
                          ( avgPosCost(i-1)*abs(prevPos(i))+ tradeda.ShortPrice(i)*abs(tradeda.ShortQty(i)))/abs(curPos(i));
                   case 'decrease'
                       avgPosCost(i) = avgPosCost(i-1);
                   case 'flatten'
                       avgPosCost(i) = avgPosCost(i-1);
               end
            end
            
            % 计算 realizedPL
            temp_idx = find(action=='decrease'|action=='flatten');
            realizedPL(temp_idx) = ( avgPosCost(temp_idx-1) - tradeda.ShortPrice(temp_idx) ).* abs(tradeda.ShortQty(temp_idx))...
                *portfolio.Instruments.(char(symbol)).Multiplier ;
            netPL = realizedPL - tradeda.TxnFee ;
            res = [tradeda, table(action,prevPos,curPos,avgPosCost,netPL,realizedPL)] ;
        end
        
        function  calcPerTradeStats(portfolio)
            
            symbols = fieldnames(portfolio.Instruments);
            nsymbols = length(symbols);
            
            for i=1:nsymbols
               
                portfolio.calcPerTradeStats_s(symbols(i));
            end
            disp(['完成了 ',num2str(nsymbols),'个对象的逐笔交易统计'])
        end
        
        function [res ] = calcPerTradeStats_s(portfolio, symbol)
           
            % 分别计算多头空头的交易
            res_long = portfolio.getPerTradeStats_long( symbol );
            res_short = portfolio.getPerTradeStats_short( symbol );
            res = [res_long;res_short];
            res = sortrows(res,'startdate','ascend');
            res.startdate = datestr(res.startdate,'yyyy-mm-dd');
            res.enddate = datestr( res.enddate,'yyyy-mm-dd');
            portfolio.perTradeStats.(char(symbol)) = res;
        end
        
        function calcAllSymbolsPL(portfolio)
            
            symbols = fieldnames(portfolio.Instruments);
            n = length(symbols);
            
            for i = 1:n
               portfolio.symbolPL.(char(symbols(i))) = portfolio.getSymbolPL(symbols(i)); 
            end
            
            disp(['完成了 ' ,num2str(n),' 个对象的交易盈亏统计' ])
        end
        
        % 将SymbolPL里各个symbol的PL数据整合成一张表
        function aggregateSymbolPL(portfolio)
            symbols = fieldnames(portfolio.symbolPL);
            n = length(symbols);
            pl = table();
            for i =1:n
               symbol = symbols(i);
               pl = [ pl; [table( symbol ),portfolio.symbolPL.(char(symbol))] ];
                
            end
            portfolio.portfolioPL.allSymbols = pl;
            
            disp('完成了各交易对象的PnL数据整合，在 portfolio.portfolioPL.allSymbols 的table表中..')
        end
        
        % 计算每个品种/交易对象的每日仓位盈亏
        function calcPositionsPL(portfolio)
            symbols = fieldnames( portfolio.positions);
            n = length(symbols);
            for i =1:n
               
                portfolio.positionsPL.(char(symbols(i))) = portfolio.calcPosPL_s( symbols(i));
                
            end
            disp([ '完成',num2str(n),'个品种/交易对象的每日仓位盈亏计算..' ])
            
        end
        % 计算单个 symbol 的每日仓位盈亏
        function [res] = calcPosPL_s(portfolio, symbol)
            position = portfolio.positions.(char(symbol));
            ndays = size(position,1);
            dates = position.Date;
            
            txnpl = portfolio.txnPL.(char(symbol));
            % 仓位当日的浮动盈亏即未实现收益
            posPL = zeros(ndays,1);
            % 每天实现的盈亏
            realizedPL = zeros(ndays,1);
            % 每天的净盈亏
            netPL = zeros(ndays,1);
            avgPosCost_long = zeros(ndays,1);
            avgPosCost_short= zeros(ndays,1);
            for i =1:ndays
                date = dates(i);
                
                %当天的交易
                txns = txnpl( txnpl.Date == date,: );
                
                % 当天的持仓
                pos = position( position.Date ==date,: );
                % 判断是否有交易数据
                if isempty(txns)
                   hastxns = 0;
                else
                   hastxns =1;
                   txns_long = txns( txns.LongQty ~= 0,: );
                   txns_short = txns( txns.ShortQty ~= 0,: );
                end
                
                if hastxns
                    if isempty(txns_long)
                        if i==1
                        avgPosCost_long(i) = 0;
                        else
                        avgPosCost_long(i) = avgPosCost_long(i-1);
                        end
                    else
                        avgPosCost_long(i) = txns_long.avgPosCost(end);    
                    end
                   if isempty(txns_short)
                       if i==1
                        avgPosCost_short(i) = 0;
                       else
                        avgPosCost_short(i) = avgPosCost_short(i-1);
                       end
                   else
                        avgPosCost_short(i) = txns_short.avgPosCost(end);    
                   end
                   
                   
                   realizedPL(i) = sum(txns.realizedPL);
                   
                   netPL(i) = sum(txns.netPL);
                   posPL(i) = ( (pos.StlPrice(1)-avgPosCost_long(i))*pos.LongPos(1)+...
                       ( avgPosCost_short(i) - pos.StlPrice(1) )*abs( pos.ShortPos(1) ) )*...
                       portfolio.Instruments.(char(symbol)).Multiplier;
                else
                   avgPosCost_long(i)=avgPosCost_long(i-1);
                   avgPosCost_short(i)=avgPosCost_short(i-1);
                   posPL(i) = ( (pos.StlPrice(1)-avgPosCost_long(i))*pos.LongPos(1)+...
                       ( avgPosCost_short(i) - pos.StlPrice(1) )*abs( pos.ShortPos(1) ) )*...
                       portfolio.Instruments.(char(symbol)).Multiplier;
                   realizedPL(i) = 0;
                   
                   netPL(i) = 0;
                end
              
                
            end
            
            cumNetPL = cumsum(netPL)+posPL;
            cumPL = cumsum( realizedPL)+posPL ;
            res = table(posPL,netPL,realizedPL,cumPL,cumNetPL,avgPosCost_long,avgPosCost_short);
            res = [position,res];
        end
        
        %-- get functions --%
        
        %计算每个 Symbol 的总体收益，将 perTradeStats的数据加总
        function [res ] = getSymbolPL(portfolio, symbol)
            
            datatable = portfolio.perTradeStats.(char(symbol));
            realizedPL = sum(datatable.realizedPL);
            netPL = sum(datatable.netPL);
            txnFee = sum(datatable.txnFee);
            nTxns = sum(datatable.nTxns);
            nTrades = size(datatable,1);
            % 根据交易盈亏计算胜率
            win = find(datatable.netPL>0);
            loss = find(datatable.netPL<=0);
            winRate = length(win)/(length(loss)+length(win));
            avgGain = mean(datatable.netPL(win));
            avgLoss = mean( datatable.netPL(loss));
            if avgLoss ~= 0 && ~isnan(avgLoss)
                plRatio = -avgGain/avgLoss;
            
            else
               plRatio = NaN; 
            end
            res = struct('realizedPL',realizedPL,'netPL',netPL,'txnFee',txnFee,'nTrades',nTrades,'nTxns',nTxns,'winrate',winRate,...
                'avgGain',avgGain,'avgLoss',avgLoss,'plRatio',plRatio);
            res = struct2table(res);
        end
        
        function [res ] = getPerTradeStats_long(portfolio, symbol)
           datatable = portfolio.txnPL.(char(symbol));
           datatable = datatable( datatable.LongQty~=0,:);
           
           if size(datatable,1) == 0
               res = table();
               return
           end
           openidx = find(datatable.action=='open');
           flatidx = find(datatable.action=='flatten'); 
           % 检查是否有 open trade,保证 opendix flatidx长度相等
           % hasopentrade = 0;
           if length(openidx) ~= length(flatidx)
               if length(openidx) - length(flatidx) == 1
                   % hasopentrade = 1;
                   warning('has open trade')
                   openidx = openidx(1:end-1);
               end
           end
           
           ntrades = length(openidx);
           startdate = zeros( ntrades,1);
           enddate = zeros( ntrades,1);
           initPos = zeros( ntrades,1);
           netPL = zeros( ntrades,1);
           realizedPL = zeros( ntrades,1);
           txnFee = zeros( ntrades,1);
           nTxns = zeros( ntrades,1);
           maxPos = zeros( ntrades,1);
           for i=1:ntrades
              startdate(i) = datatable.Date(openidx(i));
              enddate(i) = datatable.Date(flatidx(i));
              initPos(i) = datatable.curPos(openidx(i));
              maxPos(i) = max(abs(datatable.curPos( openidx(i):flatidx(i) ))) ;
              netPL(i) = sum( datatable.netPL( openidx(i):flatidx(i) ) );
              realizedPL(i) = sum( datatable.realizedPL( openidx(i):flatidx(i) ) );
              txnFee(i) = sum( datatable.TxnFee( openidx(i):flatidx(i) ) );
              nTxns(i) = flatidx(i)-openidx(i)+1;
           end
           res = table( startdate, enddate, initPos, maxPos, netPL, realizedPL, txnFee, nTxns );
           
           
           
        end
        
        function [res ] = getPerTradeStats_short(portfolio, symbol)
           datatable = portfolio.txnPL.(char(symbol));
           datatable = datatable( datatable.ShortQty~=0,:);
           if size(datatable,1) == 0
               res = table();
               return
           end
           openidx = find(datatable.action=='open');
           flatidx = find(datatable.action=='flatten'); 
           % 检查是否有 open trade,保证 opendix flatidx长度相等
           % hasopentrade = 0;
           if length(openidx) ~= length(flatidx)
               if length(openidx) - length(flatidx) == 1
                   % hasopentrade = 1;
                   warning('has open trade')
                   openidx = openidx(1:end-1);
               end
           end
           
           ntrades = length(openidx);
           startdate = zeros( ntrades,1);
           enddate = zeros( ntrades,1);
           initPos = zeros( ntrades,1);
           netPL = zeros( ntrades,1);
           realizedPL = zeros( ntrades,1);
           txnFee = zeros( ntrades,1);
           nTxns = zeros( ntrades,1);
           maxPos = zeros( ntrades,1);
           for i=1:ntrades
              startdate(i) = datatable.Date(openidx(i));
              enddate(i) = datatable.Date(flatidx(i));
              initPos(i) = datatable.curPos(openidx(i));
              maxPos(i) = min(datatable.curPos( openidx(i):flatidx(i) )) ;
              netPL(i) = sum( datatable.netPL( openidx(i):flatidx(i) ) );
              realizedPL(i) = sum( datatable.realizedPL( openidx(i):flatidx(i) ) );
              txnFee(i) = sum( datatable.TxnFee( openidx(i):flatidx(i) ) );
              nTxns(i) = flatidx(i)-openidx(i)+1;
           end
           res = table( startdate, enddate, initPos, maxPos, netPL, realizedPL, txnFee, nTxns );
            
        end
       
        
        % 计算平均年化收益率
        function [ annRtn ] = getAnnReturn(portfolio, period)
            % 默认值 period default to be 'monthly'
            if nargin == 1
               period = 'monthly'; 
            end
            if ~ismember(period,{'daily','weekly','monthly','yearly'})
                error('period的输入值只能是 daily,weekly,monthly,yearly')
            end
            
            rtn = fts2mat(portfolio.returns.(char(period)));
            switch period
                case 'daily'
                    m = mean(rtn);
                    annRtn = (1+m)^250-1;
                case 'monthly'
                    m = mean(rtn);
                    annRtn = (1+m)^12 -1;
                case 'weekly'
                    m = mean(rtn);
                    annRtn = (1+m)^52 -1;
            end % end of getAnnReturn
        end
        
        function [upstd] = getUpsideStd(portfolio, MAR ,period)
            if nargin ==1
               MAR = 0;
               period = 'daily';
            end
            if nargin == 2
               period = 'daily'; 
            end
            if ~ismember(period,{'daily','weekly','monthly','yearly'})
                error('period的输入值只能是 daily,weekly,monthly,yearly')
            end
            
            % 调整目标收益率
            T = 250;
            switch period
                case 'daily'
                    T = 250;
                    
                case 'weekly'
                    T = 52;
                    
                case 'monthly'
                    T = 12;
                    
            end
            target = ( 1+MAR )^( 1/T ) - 1 ;
            
            rtn = fts2mat(portfolio.returns.(char(period)));
            rtn_posi = rtn( rtn>target );
            upstd = std(rtn_posi);
        end % end of getUpsideStd
        
        function [downstd] = getDownsideStd(portfolio,MAR, period)
            if nargin ==1
               MAR = 0;
               period = 'daily';
            end
            if nargin == 2
               period = 'daily'; 
            end
            if ~ismember(period,{'daily','weekly','monthly','yearly'})
                error('period的输入值只能是 daily,weekly,monthly,yearly')
            end
            
            % 调整目标收益率
            T = 250;
            switch period
                case 'daily'
                    T = 250;
                    
                case 'weekly'
                    T = 52;
                    
                case 'monthly'
                    T = 12;
                    
            end
            target = ( 1+MAR )^( 1/T ) - 1 ;
            
            rtn = fts2mat(portfolio.returns.(char(period)));
            rtn_nega = rtn( rtn<target );
            downstd = std(rtn_nega);
        end % end of getDownsideStd
        
        function [res] = getSharpe_period(portfolio,startdate,enddate,riskfree)
            if nargin == 3
               riskfree =0; 
            end
            
            if ischar(startdate)
               startdate = datenum(startdate) ;
            end
            if ischar(enddate)
               enddate = datenum(enddate) ;
            end
            % 将risk free rate 转换为日周期值
            if riskfree ~= 0 
               riskfree = (riskfree+1)^(1/250)-1; 
            end
            
            rtn = portfolio.returns.daily(portfolio.returns.daily.dates>=startdate & portfolio.returns.daily.dates<=enddate);
            res = sharpe(fts2mat(rtn),riskfree);
        end
        
        function [res] = getSharpe_eachmonth(portfolio,riskfree)
            if nargin == 1
               riskfree = 0; 
            end
            % get end dates of each month. a datenum vector
            enddates = portfolio.returns.monthly.dates;
            startdates = enddates;
            nmon = length(enddates);
            v = datevec(startdates);
            v(:,3) = 1;
            startdates = datenum(v);
            
            shp = zeros(nmon,1);
            for i =1:nmon
                shp(i) = portfolio.getSharpe_period(startdates(i),enddates(i),riskfree);
            end
            res = fints( enddates,shp, 'sharpe_eachmonth','monthly') ; 
            
        end
        %-- end of get functions --%
        
        %-- %
        
        
        
        %-- %
        
        %- chart functions -%
        
       
        
        % 画收益率分布图
        function chartReturnDistribution(portfolio,freq)
            % period的默认值是daily
            if nargin == 1
               freq = 'daily'; 
            end
            if ~ismember(freq,{'daily','weekly','monthly'})
                error('period的输入值只能是 daily,weekly,monthly')
            end
            n = 50;
            if length( portfolio.returns.(char(freq)) ) <250
               n = round( length( portfolio.returns.(char(freq)) )/5 );
            end
            figure
            histfit(fts2mat(portfolio.returns.(char(freq))),n)
            switch freq
                case 'daily'
                    d = '日';
                case 'weekly'
                    d = '周';
                case 'monthly'
                    d = '月';
            end
            title([d,'收益率分布'])
        end
        % 画月度收益的直方图
        % monthly return bar chart
        function chartMonthlyReturn(portfolio)
            figure
            bar(portfolio.returns.monthly.dates,fts2mat(portfolio.returns.monthly))
            title('月度收益率')
            xlabel('月/年')
            ylabel('收益率')
            dateaxis('x',12)
        end
        
        function chartPosition(portfolio, symbol)
           pospl = portfolio.positionsPL.(char(symbol));
           dates = pospl.Date;
           position_long = pospl.LongPos;
           position_short = pospl.ShortPos;
           cumPL = pospl.cumPL;
           cumNetPL = pospl.cumNetPL;
           
           figure
           subplot(3,1,1)
           plot(dates,pospl.StlPrice)
           title([symbol,'的Settle Price'])
           legend('Settle price')
           dateaxis('x',2);
           
           subplot(3,1,2)
           plot(dates, cumPL)
           hold on 
           plot(dates, cumNetPL)
           legend('cumPL','cumPL.net')
           title([symbol,'的PnL'])
           dateaxis('x',2)
           hold off
           
           subplot(3,1,3)
           stairs(dates,position_long)
           hold on
           stairs(dates,position_short)
           title('position')
           dateaxis('x',2)
           hold off
           
           
        end
        
        function chartSharpe_eachmonth(portfolio,riskfree)
            if nargin == 1
               riskfree = 0; 
            end
           shps = portfolio.getSharpe_eachmonth(riskfree);
            figure
            bar(shps.dates,fts2mat(shps))
            title('每个月的日收益率的sharpe ratio')
            dateaxis('x',12)
        end
        
        function chartEquity(portfolio)
           eqty = portfolio.equity;
            figure
            plot(eqty.dates,fts2mat(eqty))
            dateaxis('x',2)
            legend('Equity')
            title('每日权益曲线')
            
        end
        
        function chartDrawdown(portfolio)
            eqty = fts2mat(portfolio.equity);
            dates = portfolio.equity.dates;
            cmmax = cummax(eqty);
            dd = cmmax - eqty;
            ddrate = dd./cmmax ;
            
            figure
            subplot(2,1,1)
            plot(dates,eqty)
            dateaxis('x',2)
            legend('Equity')
            title('每日权益曲线')
            subplot(2,1,2)
            plotyy(dates,dd,dates,ddrate)
            legend('drawdown','drawdown rate')
            title('回撤')
            dateaxis('x',2)
            
        end
        %-- end of chart functions --%
        
        %-- disp functions
        function dispAllSymbolPL(portfolio)
           datatable = portfolio.portfolioPL.allSymbols;
           datatable = sortrows(datatable,'netPL','descend');
           disp(datatable)
        end
        
        
        
        %-- end of disp functions
        
        %-- output functions 
        function outAllSymbolPL(portfolio)
           datatable = portfolio.portfolioPL.allSymbols;
           datatable = sortrows(datatable,'netPL','descend');
           
           writetable(datatable,'allSymbolPnL.csv','Delimiter',',')
           disp('完成输出')
           
        end
        
        function outBasicStats(portfolio)
            fields = fieldnames(portfolio.basicStats);
            n = length(fields);
            temp_table = table();
            for i=1:n
                
                temp_table = [temp_table; portfolio.basicStats.( char( fields(i)) ) ] ;
            end
            period = fields;
            temp_table = [ table(period), temp_table];
            writetable(temp_table,'basicStats.csv','Delimiter',',')
            disp('完成输出')
        end
        
        function outRatios(portfolio)
           
            t1 = struct2table(portfolio.performance.daily);
            t2 = struct2table(portfolio.performance.weekly);
            t3 = struct2table(portfolio.performance.monthly);
            names = { 'sharpe','sortino','sterling','sterling_maxDD','burke','burke_maxDD','upsidePotential','sharpe_upside','sharpe_downside' };
            t1 = t1(:,names);
            t2 = t2(:,names);
            t3 = t3(:,names);
            tab = [t1;t2;t3];
            
            writetable(tab,'ratios.csv','Delimiter',',')
            disp('完成输出')
        end
        
        % 输出所有的开仓-平仓交易统计数据
        function outPerTradeStats(portfolio)
           
            symbols = fieldnames(portfolio.perTradeStats);
            n = length(symbols);
            tab = table();
            for i=1:n
               temp = portfolio.perTradeStats.(char(symbols(i)));
               s = size(temp,1);
               symbol = zeros(s,1);
               symbol = cellstr(char(symbol));
               for j = 1:s
                  symbol(j)=symbols(i); 
               end
               temp = [ table( symbol ),temp];
               tab = [tab;temp];
            end
            writetable(tab,'perTradeStats.csv','Delimiter',',')
            disp('完成输出')
        end
        
        function outDrawdownStats(portfolio)
           datatable =  portfolio.drawdownStats;
            writetable(datatable,'drawdownStats.csv','Delimiter',',')
            disp('完成输出')
        end
        
        %-- end of output functions
        
        % 检查是否同时在同一个品种上持有多空仓
        function checkLongShortPos( portfolio )
            
            x = portfolio.PosData.LongPos ~= 0 & portfolio.PosData.ShortPos ~= 0;
            if any(x)
               warning('有品种同时持有多空仓') 
               disp(portfolio.PosData(x,:))
            end
            
        end
        
        
    end
    
end

