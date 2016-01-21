classdef mPortfolio < handle
    % mPORTFOLIO Portfolio(Ͷ�����)��
    %   Detailed explanation goes here
    
    % define properties
    % ��������
    % equity: Ȩ���ʱ������
    % positions: ��λ���ݱ��������н��׶���Ĳ�λ��Ϣ
    % trades: ��ʽ������ݱ��������н��׶���Ľ�����Ϣ
    properties
        equity
        returns
        PosData % ���н��׵Ľ��ڹ��ߵĲ�λ��Ϣ����
        TradeData % ���н��׵Ľ��ڹ��ߵĽ��׳ɽ�����
        Instruments %
        positions % ��λ���ݣ�ÿ�����ڹ��߷��䵥����һ��table
        trades % �ɽ����ݣ�ÿ�����ڹ��߷��䵥����һ��table
        drawdownStats % �س�����
        perTradeStats % ����ӿ��ֵ�ƽ�ֵ���ʽ�������ͳ��
        txnPL % ����ÿ�ʳɽ����ݵ�ӯ������
        positionsPL %����ÿ�� symbol �Ĳ�λӯ��������
        symbolPL % ����ÿ�� Symbol �ϵĺϼ�ӯ��
        portfolioPL % ����Ͷ����ϲ����ӯ������
        performance % ����ҵ�������ĸ��� ratio
        basicStats % �����ʵĻ���ͳ����
        
    end
    
    methods
        % constructor
        function portf = mPortfolio(equity,posData,tradeData)
            % ��������posData��tradeData�е�symbol�Ƿ�һ��
            % ��¼��λ�ͳɽ���¼������һ���Եĵ�һ���Ǽ���������Ľ��ڹ��ߵ����֣���Symbol�ֶΣ��Ƿ���ȫһ��
            % ��λ�����ǿ���ͨ���ɽ���ϸ������ģ����Դӳɽ���ϸ����ó��Ĳ�λ���ݿ��Ժ�posData�Ƚ��Լ��posData��tradeData��һ����
            % �˴�ֻ����Symbol�ֶεļ��
            % ���tradeData��posData�е�symbol������һ�£�˵�����ݿ����д�©������tradeData�а��������ڽ�������
            if length( unique( posData.Symbol ) ) ~= length( unique( tradeData.Symbol ) )
            warning('posData��tradeData�е�symbol������һ��,����tradeData�������ڽ�������,��posDataΪÿ�ս����Ĳ�λ')
            end
            
            % ���Եĳ�ʼ��
            portf.equity = equity; %ÿ��Ȩ��
            portf.PosData = posData;
            portf.TradeData = tradeData;
            
            portf.Instruments = struct(); % ���ڴ���ÿ��instrument�������� Multiplier
            portf.positions = struct(); % ���ڴ���ÿ��instrument��position data
            portf.trades = struct(); % ���ڴ���ÿ��instrument��trade data
            portf.performance = struct(); % ���ڴ���ҵ�������Ľ��
            portf.basicStats = struct(); % ���ڴ����������еĻ���ͳ����
            portf.returns = struct();
            portf.drawdownStats = table();
            portf.perTradeStats = struct();
            portf.txnPL = struct();
            portf.symbolPL = struct();
            portf.portfolioPL = struct();
            portf.positionsPL = struct();
            
            disp( '�ɹ�����mPortfolio�����' )
        end
        
        % ��portfolio��������׵Ľ��ڹ��ߵ���Ϣ���ݣ���tradeData������
        function addInstruments(portfolio, varargin)
            if isempty(portfolio.TradeData)
                error('TradeData is empty, cannot create instruments from TradeData. TradeData�ǿյ�')
            end
            if nargin == 1
                error('Please input the variable we need in TradeData. �������ֶ���') 
            end
            
            labels = { 'Symbol',char(varargin)};
            [Lia,~] = ismember( labels,fieldnames(portfolio.TradeData) );
            if ~all(Lia)   
                error('�����ֵ���еĲ���TradeData���ֶ�')
            end
            temp = unique(portfolio.TradeData(:,labels) );
            symbolnames = unique(temp.Symbol);
            if size(temp,1) ~= length( symbolnames )
                error('���ڽ��ڹ����в�ͬ���������ݣ�����TradeData')
            end
            for i=1:length(symbolnames)
               
               portfolio.addInstrument(temp.Symbol(i),temp(i,varargin) );
            end
            
            disp( ['������',num2str(length(symbolnames)),'��Instrument����'] )
        end % end of function addInstruments
        
        
        % ��portfolio��������׵Ľ��ڹ��ߵ���Ϣ���ݣ������Զ�����������
        function addInstrument(portfolio, symbol, da)
            symbol = char(symbol);
            portfolio.Instruments.(symbol) = da;
        end
        
        % create positions
        % ���� positions �еĲ�λ���ݱ�ÿ��instrument����һ�������ı�
        % �� PosData �а�symbol�ָ����ݱ�
        function cpositions(portfolio)
            symbolnames = unique( portfolio.PosData.Symbol );
            for i = 1:length(symbolnames)
               symbol = symbolnames(i);
              %idx = find(portfolio.PosData.Symbol==symbol);
               portfolio.positions.(char(symbol)) = portfolio.PosData(portfolio.PosData.Symbol==symbol,:);
            end
            
            disp( ['������',num2str(length(symbolnames)),'��position���ݱ�'] )
        end % end of cpositions
        
        % create trades
        % ���� trades �еĳɽ���¼���ݱ�ÿ��instrument����һ�������ı�
        % �� TradeData �а�symbol�ָ����ݱ�
        function ctrades(portfolio)
            symbolnames = unique( portfolio.TradeData.Symbol );
            for i = 1:length(symbolnames)
               symbol = symbolnames(i);
              %idx = find(portfolio.PosData.Symbol==symbol);
               portfolio.trades.(char(symbol)) = portfolio.TradeData(portfolio.TradeData.Symbol==symbol,:);
            end
            
            disp( ['������',num2str(length(symbolnames)),'��trade���ݱ�'] )
        end % end of ctrades
        
         % ����ͳ����
        function calcBasicStats(portfolio,period)
            % period��Ĭ��ֵ��daily
            if nargin == 1
               period = 'daily'; 
            end
            if ~ismember(period,{'daily','weekly','monthly'})
                error('period������ֵֻ���� daily,weekly,monthly')
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
            disp( [period,' return��basicStats��'])
            disp(portfolio.basicStats.(period))
            
        end % end of calcBasicStats
        
        % ����������ʣ�����equity�����е�equity���е�Ȩ��ֵ���㣬���ڵ���Ȩ��ֵ���е�����
        % returns���еĵ�һ��ֵΪ NaN
        function calcDailyReturn(portfolio,method)
            % method��Ĭ��ֵ��Simple
            if nargin == 1
               method = 'Simple'; 
            end
             % �����Ȩ�������Ƿ�Ϊ��
            if isempty(portfolio.equity)
               error('equity is empty. equityΪ��.') 
            end
            % ��� method ��ֵ�Ƿ�ΪSimple��Continuous
            if ~ismember(method,{'Simple','Continuous'})
                error('method������ֵֻ���� Simple �� Continuous')
            end
             portfolio.returns.daily = tick2ret(portfolio.equity.equity,'method',method) ;
             
        end
        
        function calcWeeklyReturn(portfolio,method)
            % method��Ĭ��ֵ��Simple
            if nargin == 1
               method = 'Simple'; 
            end
            % �����Ȩ�������Ƿ�Ϊ��
            if isempty(portfolio.equity)
               error('equity is empty. equityΪ��.') 
            end
            % ��� method ��ֵ�Ƿ�ΪSimple��Continuous
            if ~ismember(method,{'Simple','Continuous'})
                error('method������ֵֻ���� Simple �� Continuous')
            end
            % ����Ȩ������ת��Ϊ�ܶ�����
            equity_weekly = toweekly(portfolio.equity);
            % ������������
            portfolio.returns.weekly = tick2ret(equity_weekly,'method',method);
        end
        
        function calcMonthlyReturn(portfolio,method)
            % method��Ĭ��ֵ��Simple
            if nargin == 1
               method = 'Simple'; 
            end
            % �����Ȩ�������Ƿ�Ϊ��
            if isempty(portfolio.equity)
               error('equity is empty. equityΪ��.') 
            end
            % ��� method ��ֵ�Ƿ�ΪSimple��Continuous
            if ~ismember(method,{'Simple','Continuous'})
                error('method������ֵֻ���� Simple �� Continuous')
            end
            % ����Ȩ������ת��Ϊ�¶�����
            equity_monthly = tomonthly(portfolio.equity);
            % ������������
            portfolio.returns.monthly = tick2ret(equity_monthly,'method',method);
        end % end of calcMonthlyReturn
        
        % �������ձ��� sharpe ratio
        % riskfree Ϊ�޷�������, default 0
        % period Ϊ�����ʵ�����, default 'daily'
        function [shpr] = calcSharpe(portfolio,riskfree,period)
            if nargin == 1
               riskfree = 0; % risk free rate default to be zero
               period = 'daily'; % period default to be daily 
            end
            if nargin == 2
               period = 'daily'; % period default to be daily 
            end
            if ~ismember(period,{'daily','weekly','monthly','yearly'})
                error('period������ֵֻ���� daily,weekly,monthly,yearly')
            end
            
            % �� risk free rate ת����ʹ���ں� return ������һ��
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
            
            % �����������ѡ�������ʽ�������ú���ʱ�������������������������ֵ�������ڿ���̨��ӡ������
            if nargout == 1
            shpr = portfolio.performance.(char(period)).sharpe;
            else 
            disp( [ period,' return��sharpe ratioΪ', num2str(portfolio.performance.(char(period)).sharpe) ] )
            end
        end
        
        function calcMaxDrawdown(portfolio)
            portfolio.performance.daily.maxDrawdown = maxdrawdown(fts2mat(portfolio.equity.equity),'return');
            disp( [ 'MaxDrawdownΪ', num2str(portfolio.performance.daily.maxDrawdown) ] )
        end
        
        % ���� sortino ratio
        % MAR: �껯��Ŀ�������ʣ�Ͷ���߳�������Ϊ 0, ��ʷ��ֵ������, �޷�������, ���׼����
        % �����µ���ƫ��� downsideDeviation(x,target)
        % target Ϊ��MAR����Ϊ��x�����ڶ�Ӧ��ֵ. ���磺x ���������ʣ���targetΪ���������Ŀ��������
        function [sortino] = calcSortino(portfolio, MAR, period)
            % set default values, ����Ĭ��ֵ
            if nargin == 1
               MAR = 0;
               period = 'monthly';
            end
            if nargin == 2
               period = 'monthly'; 
            end
            % ����Ŀ��������
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
                disp( [ period,' return����',num2str(MAR),'ΪĿ�������ʵ�Sortino ratioΪ', num2str(portfolio.performance.(char(period)).sortino) ] )
            end
        end % end of calcSortino
        
        % upside potential ratio
        % ��ЩͶ����Ѱ�ҵĲ���һ��ʱ����ƽ��������ߵľ���������Щ����MAR��ƽ����������ߵľ���
        % Sortino, van der Meer and Plantinga(1999)
        % ����������пռ�������Sortino�����еķ����еĶ�������
        % ���пռ�ָ���ǳ���MAR��Ԥ�������ʣ����Կ����ǳɹ���Ǳ��
        % upside potential ratio = E[ r - MAR ]/downside risk
        % �� Sortino һ����upside potential ratio Ҳ���껯������
        function [uppr] = calcUpsidePotential(portfolio,MAR,period)
             % set default values, ����Ĭ��ֵ
            if nargin == 1
               MAR = 0;
               period = 'monthly';
            end
            if nargin == 2
               period = 'monthly'; 
            end
             % ����Ŀ��������
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
                disp( [ 'Upside potential ratio ���� ',period,' return�� ',num2str(uppr) ] )
            end
        end
        
        
        % Sterling ratio ����Sterling����
        % Sterling��Burke�����ܵ�CTAs�Ĺ㷺��ӭ����Ϊ���������ʳ���������������������õ��Ƿ��棺
        % ��������󻯣��������ϸ������ʧ.
        % Sterling���ʿ����ûس����������գ����Ա�Sortino���ʸ�����һ��
        % sterling = (mean(r) - riskfree )/mean(drawdown)
        % ��ĸ������Ϊ��ֵ������һЩ�س���ƽ��ֵ��"��ֵ����"���д�����.�˴����䶨��Ϊ���лس���ƽ��ֵ
        function [sterling] = calcSterling(portfolio,riskfree,period)
            
            if nargin == 1
               riskfree = 0; % risk free rate default to be zero
               period = 'monthly'; % period default to be daily 
            end
            if nargin == 2
               period = 'monthly'; % period default to be daily 
            end
            if ~ismember(period,{'daily','weekly','monthly','yearly'})
                error('period������ֵֻ���� daily,weekly,monthly,yearly')
            end
         
            eqty = fts2mat(portfolio.equity);
            eqty = eqty(:,1);
            % cumulative maximum
            % index of cmmax is the same as eqty
            cmmax = cummax(eqty,1); 
            % �س����㹫ʽ�� drawdown :=  cummax - eqty (which is positive)
            % dd and ddrate have same indexing as eqty
            dd =  cmmax - eqty ;
            ddrate = dd ./cmmax ;
            
            
            annRtn = portfolio.getAnnReturn( period );
            mdwn = mean( ddrate );
            sterling = ( annRtn - riskfree)/mdwn;
            portfolio.performance.(char(period)).sterling = sterling;
        end % end of calcSterling
        
        
        %�ø����س��ڵ����س��������еĻس�
        function [sterling_maxDD] = calcSterling_maxDD(portfolio,riskfree,period)
            
            if nargin == 1
               riskfree = 0; % risk free rate default to be zero
               period = 'monthly'; % period default to be daily 
            end
            if nargin == 2
               period = 'monthly'; % period default to be daily 
            end
            if ~ismember(period,{'daily','weekly','monthly','yearly'})
                error('period������ֵֻ���� daily,weekly,monthly,yearly')
            end
  
            annRtn = portfolio.getAnnReturn(period);
            mdwn = mean( portfolio.drawdownStats.maxDrawdownPct./100 );
            sterling_maxDD = ( annRtn - riskfree)/mdwn;
            portfolio.performance.(char(period)).sterling_maxDD = sterling_maxDD;
        end % end of calcSterling_maxDD
        
        % Burke ratio
        % Burke������Sterling ratio,
        % Burke(1994)����ʹ��ÿ���س���ƽ���͵�ƽ�����Գͷ�������ºͻس���Եķ��Ƚϴ�Ļس�
        function [ burke ] = calcBurke(portfolio, riskfree,period)
            if nargin == 1
               riskfree = 0; % risk free rate default to be zero
               period = 'monthly'; % period default to be daily 
            end
            if nargin == 2
               period = 'monthly'; % period default to be daily 
            end
            if ~ismember(period,{'daily','weekly','monthly','yearly'})
                error('period������ֵֻ���� daily,weekly,monthly,yearly')
            end
          
            eqty = fts2mat(portfolio.equity);
            eqty = eqty(:,1);
            % cumulative maximum
            % index of cmmax is the same as eqty
            cmmax = cummax(eqty,1); 
            % �س����㹫ʽ�� drawdown :=  cummax - eqty
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
                error('period������ֵֻ���� daily,weekly,monthly,yearly')
            end
         
            dwn = portfolio.drawdownStats.maxDrawdownPct./100 ;
            deno = sqrt( sum(dwn.^2) );
            burke_maxDD = (portfolio.getAnnReturn(period) - riskfree)/deno;
            portfolio.performance.(char(period)).burke_maxDD = burke_maxDD;
        end % end of calcBurke_maxDD
        
        
        
        % ͳ��ÿһ���س��ڵ�����
        function calcDrawdownStats(portfolio)
            % �����Ȩ�������Ƿ�Ϊ��
            if isempty(portfolio.equity)
               error('equity is empty. equityΪ��.') 
            end
            % convert fts to array
            eqty = fts2mat(portfolio.equity);
            eqty = eqty(:,1);
            % cumulative maximum
            % index of cmmax is the same as eqty
            cmmax = cummax(eqty,1); 
            % �س����㹫ʽ�� drawdown :=  cummax - eqty
            % dd and ddrate have same indexing as eqty
            dd =  cmmax - eqty ;
            ddrate = dd ./cmmax ;
            % �س��յ�index. get the index number of drawdown days
            % for indexing in eqty
            ddidx = find(dd>0);
          
            startidx = [ddidx(1);ddidx( [1;diff(ddidx)]>1 )];
            endidx = [ddidx( [diff(ddidx);1]>1 ) ; ddidx(end)] ;
            startdate = cellstr(datestr(portfolio.equity.dates(startidx),'yyyy-mm-dd'));
            enddate =  cellstr(datestr(portfolio.equity.dates(endidx),'yyyy-mm-dd'));
            timespan = endidx - startidx +1;
            %-% ÿ���س����ڿ����ж���ֲ��ߵ㣬Ѱ��local maximum�ķ���������
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
            disp('����س��������,�ɵ���drawdownStats���Բ鿴....')
        end
        
        % �����µ�Ƶ��
        function [dwnfreq] = calcDownsideFreq(portfolio,MAR,period)
            if nargin ==1
               MAR = 0;
               period = 'monthly';
            end
            if nargin == 2
               period = 'monthly'; 
            end
            if ~ismember(period,{'daily','weekly','monthly','yearly'})
                error('period������ֵֻ���� daily,weekly,monthly,yearly')
            end
            
             % ����Ŀ��������
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
        
        % ���� �������sharpe
        function [sharpe_upside] = calcSharpe_upside(portfolio,MAR,period)
            
              % set default values, ����Ĭ��ֵ
            if nargin == 1
               MAR = 0;
               period = 'monthly';
            end
            if nargin == 2
               period = 'monthly'; 
            end
             % ����Ŀ��������
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
                          % set default values, ����Ĭ��ֵ
            if nargin == 1
               MAR = 0;
               period = 'monthly';
            end
            if nargin == 2
               period = 'monthly'; 
            end
             % ����Ŀ��������
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
            disp(['����� ',num2str(nsymbols),'������Ľ������ݵļ���'])
        end
        
        function [ result ] = calcTxn_s(portfolio,symbol)
            n = size(portfolio.trades.(char(symbol)),1);
            portfolio.trades.(char(symbol)).txnId = (1:n)' ;
            % �����ͷ
            txn_long = portfolio.calcTxn_long(symbol);
            
            % �����ͷ
            txn_short = portfolio.calcTxn_short(symbol);
            % �ϲ����
            result = [txn_long;txn_short];
            result = sortrows(result,'txnId','ascend');
            % ���� ��ճֲ�
            longPos = cumsum(result.LongQty);
            longPos = table(longPos);
            shortPos = cumsum(result.ShortQty);
            shortPos = table(shortPos);
            result = [result,longPos,shortPos];
            % ���� cumPL
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
            
            % ���� current position and previous position
            curPos = cumsum(tradeda.LongQty);
            % previous position of first record is default to be 0, unless
            % we need to consider situation that the portfolio has position
            % before the performence analysis period. and this will be
            % continued to develop in the future.
            prevPos(1) = 0;
            prevPos(2:end) = curPos(1:end-1);
            
            % ���� action
            action(prevPos==0)='open';
            action( curPos==0 ) ='flatten';
            action( prevPos~=0 & curPos ~=0 & curPos > prevPos)='increase';
            action( prevPos~=0 & curPos ~=0 & curPos < prevPos)='decrease';
            
            % ����ÿ�ʽ���֮���ƽ���ֲֳɱ� 
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
            
            % ���� realizedPL
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
            
            % ���� action
            action(prevPos==0) = 'open';
            action( curPos==0 ) = 'flatten';
            action( prevPos~=0 & curPos ~=0 & curPos > prevPos)='decrease';
            action( prevPos~=0 & curPos ~=0 & curPos < prevPos)='increase';
            
            % ����ÿ�ʽ���֮���ƽ���ֲֳɱ� 
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
            
            % ���� realizedPL
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
            disp(['����� ',num2str(nsymbols),'���������ʽ���ͳ��'])
        end
        
        function [res ] = calcPerTradeStats_s(portfolio, symbol)
           
            % �ֱ�����ͷ��ͷ�Ľ���
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
            
            disp(['����� ' ,num2str(n),' ������Ľ���ӯ��ͳ��' ])
        end
        
        % ��SymbolPL�����symbol��PL�������ϳ�һ�ű�
        function aggregateSymbolPL(portfolio)
            symbols = fieldnames(portfolio.symbolPL);
            n = length(symbols);
            pl = table();
            for i =1:n
               symbol = symbols(i);
               pl = [ pl; [table( symbol ),portfolio.symbolPL.(char(symbol))] ];
                
            end
            portfolio.portfolioPL.allSymbols = pl;
            
            disp('����˸����׶����PnL�������ϣ��� portfolio.portfolioPL.allSymbols ��table����..')
        end
        
        % ����ÿ��Ʒ��/���׶����ÿ�ղ�λӯ��
        function calcPositionsPL(portfolio)
            symbols = fieldnames( portfolio.positions);
            n = length(symbols);
            for i =1:n
               
                portfolio.positionsPL.(char(symbols(i))) = portfolio.calcPosPL_s( symbols(i));
                
            end
            disp([ '���',num2str(n),'��Ʒ��/���׶����ÿ�ղ�λӯ������..' ])
            
        end
        % ���㵥�� symbol ��ÿ�ղ�λӯ��
        function [res] = calcPosPL_s(portfolio, symbol)
            position = portfolio.positions.(char(symbol));
            ndays = size(position,1);
            dates = position.Date;
            
            txnpl = portfolio.txnPL.(char(symbol));
            % ��λ���յĸ���ӯ����δʵ������
            posPL = zeros(ndays,1);
            % ÿ��ʵ�ֵ�ӯ��
            realizedPL = zeros(ndays,1);
            % ÿ��ľ�ӯ��
            netPL = zeros(ndays,1);
            avgPosCost_long = zeros(ndays,1);
            avgPosCost_short= zeros(ndays,1);
            for i =1:ndays
                date = dates(i);
                
                %����Ľ���
                txns = txnpl( txnpl.Date == date,: );
                
                % ����ĳֲ�
                pos = position( position.Date ==date,: );
                % �ж��Ƿ��н�������
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
        
        %����ÿ�� Symbol ���������棬�� perTradeStats�����ݼ���
        function [res ] = getSymbolPL(portfolio, symbol)
            
            datatable = portfolio.perTradeStats.(char(symbol));
            realizedPL = sum(datatable.realizedPL);
            netPL = sum(datatable.netPL);
            txnFee = sum(datatable.txnFee);
            nTxns = sum(datatable.nTxns);
            nTrades = size(datatable,1);
            % ���ݽ���ӯ������ʤ��
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
           % ����Ƿ��� open trade,��֤ opendix flatidx�������
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
           % ����Ƿ��� open trade,��֤ opendix flatidx�������
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
       
        
        % ����ƽ���껯������
        function [ annRtn ] = getAnnReturn(portfolio, period)
            % Ĭ��ֵ period default to be 'monthly'
            if nargin == 1
               period = 'monthly'; 
            end
            if ~ismember(period,{'daily','weekly','monthly','yearly'})
                error('period������ֵֻ���� daily,weekly,monthly,yearly')
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
                error('period������ֵֻ���� daily,weekly,monthly,yearly')
            end
            
            % ����Ŀ��������
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
                error('period������ֵֻ���� daily,weekly,monthly,yearly')
            end
            
            % ����Ŀ��������
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
            % ��risk free rate ת��Ϊ������ֵ
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
        
       
        
        % �������ʷֲ�ͼ
        function chartReturnDistribution(portfolio,freq)
            % period��Ĭ��ֵ��daily
            if nargin == 1
               freq = 'daily'; 
            end
            if ~ismember(freq,{'daily','weekly','monthly'})
                error('period������ֵֻ���� daily,weekly,monthly')
            end
            n = 50;
            if length( portfolio.returns.(char(freq)) ) <250
               n = round( length( portfolio.returns.(char(freq)) )/5 );
            end
            figure
            histfit(fts2mat(portfolio.returns.(char(freq))),n)
            switch freq
                case 'daily'
                    d = '��';
                case 'weekly'
                    d = '��';
                case 'monthly'
                    d = '��';
            end
            title([d,'�����ʷֲ�'])
        end
        % ���¶������ֱ��ͼ
        % monthly return bar chart
        function chartMonthlyReturn(portfolio)
            figure
            bar(portfolio.returns.monthly.dates,fts2mat(portfolio.returns.monthly))
            title('�¶�������')
            xlabel('��/��')
            ylabel('������')
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
           title([symbol,'��Settle Price'])
           legend('Settle price')
           dateaxis('x',2);
           
           subplot(3,1,2)
           plot(dates, cumPL)
           hold on 
           plot(dates, cumNetPL)
           legend('cumPL','cumPL.net')
           title([symbol,'��PnL'])
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
            title('ÿ���µ��������ʵ�sharpe ratio')
            dateaxis('x',12)
        end
        
        function chartEquity(portfolio)
           eqty = portfolio.equity;
            figure
            plot(eqty.dates,fts2mat(eqty))
            dateaxis('x',2)
            legend('Equity')
            title('ÿ��Ȩ������')
            
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
            title('ÿ��Ȩ������')
            subplot(2,1,2)
            plotyy(dates,dd,dates,ddrate)
            legend('drawdown','drawdown rate')
            title('�س�')
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
           disp('������')
           
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
            disp('������')
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
            disp('������')
        end
        
        % ������еĿ���-ƽ�ֽ���ͳ������
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
            disp('������')
        end
        
        function outDrawdownStats(portfolio)
           datatable =  portfolio.drawdownStats;
            writetable(datatable,'drawdownStats.csv','Delimiter',',')
            disp('������')
        end
        
        %-- end of output functions
        
        % ����Ƿ�ͬʱ��ͬһ��Ʒ���ϳ��ж�ղ�
        function checkLongShortPos( portfolio )
            
            x = portfolio.PosData.LongPos ~= 0 & portfolio.PosData.ShortPos ~= 0;
            if any(x)
               warning('��Ʒ��ͬʱ���ж�ղ�') 
               disp(portfolio.PosData(x,:))
            end
            
        end
        
        
    end
    
end

