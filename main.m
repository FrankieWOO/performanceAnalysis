%% ˵��
% 1. MATLAB �汾����(�������MATLAB R2014b��)��
% 1.1 MATLAB������table�������������滻dataset����, �����´�����ʹ����table
%     �������MATLAB�汾�����޷�����table����ʹ�ù�����ͬ��dataset
% 1.2 fts (financial time series)������Ҫ Financial toolbox ��֧��
%
% 2. mPortfolio ��
%   ����������������Ƶ�˼���Զ��岢ʹ���� mPortfolio �ࣨǰ׺m���ں�Financial toolbox��portfolio������)
%   ��������Ͷ����������еļ����ͳ�Ƶĺ�����װ���� mPortfolio ���У�����mPortfolio.m�ļ��в鿴�޸�
%   ʹ�����װ���Խ����������ͻ��ʹ����ṹ������������ע��Ч�ʵ�Ӧ���У�ʱ�俪������Ҫ���ǵ�����
%   
% 3. ���д��룺 �밴section˳���������д��뵼�벢����portfolio����ҵ��������ش�����ظ�����
%% �������ݼ�
% �������ݼ�
% import datasets of position data, product data and trade data
% in recent matlab version 'table' is introduced to replace dataset,
% the function dataset is replaced by readtable.
%  
% posData = dataset('XLSFile','Position Data.xlsx');
% productData = dataset('XLSFile','Product Data.xlsx');
% tradeData = dataset('XLSFile','Position Data.xlsx');

posData = readtable('\data\Position Data.xlsx');
productData = readtable('\data\Product Data.xlsx');
tradeData = readtable('\data\Trade Data.xlsx');

%% ���õ�������ݱ�������
% ��������ɲ鿴���������ݼ���������
% these commands can access the variable names of the tables; 
% VariableNames need to be set manually.
% (dataset)
% posData.Properties.VarNames
% productData.Properties.VarNames
% tradeData.Properties.VarNames
% (table)
% 
% 
% ����������ڲ�ѯ���ݼ�����������������
% (dataset)
% datasetfun(@class,posData,'UniformOutput',false)
% (table)
% posData.Properties.VariableDescriptions

% modify the names manually
% �޸ĸ��е�����
namesPosData = {'Date','Symbol',...
    'LongPos','ShortPos','PrevStlPrice','StlPrice','Multiplier'};
posData.Properties.VariableNames([1,2,3,5,7,8,9]) = namesPosData; 
namesProductData = {'Date','Equity'};
productData.Properties.VariableNames([1,2]) = namesProductData;
namesTradeData = {'Date',...
    'Symbol','LongQty','LongPrice','ShortQty','ShortPrice','PrevStlPrice',...
    'StlPrice','Multiplier','TxnFee'};
tradeData.Properties.VariableNames([1,2,3,4,5,6,7,8,9,10]) = namesTradeData;


%% clean the Nan row in tables. ��Ե�������ݱ���п�������
% �������ݱ���DateΪ�յ���
productData = cleanNanRow(productData,'Date');
posData = cleanNanRow(posData,'Date');
tradeData = cleanNanRow(tradeData,'Date');

% ���� tradeData ��û�в�����¼����. cellstring�п�ֵ���� ''
tradeData = cleanNanRow(tradeData,'Symbol');

%% ��Ե�������ݱ������������ת��
%��datestringת��Ϊ datenum����
posData.Date = datenum(posData.Date);
productData.Date = datenum(num2str(productData.Date),'yyyymmdd');
tradeData.Date = datenum( tradeData.Date);

% ��Symbol��cell array of stringת��Ϊcategorical����
posData.Symbol = categorical(posData.Symbol);
tradeData.Symbol = categorical(tradeData.Symbol);
%% ����portfolio����

% ����Ȩ���ʱ������. create financial time series object of equity
equity = fints(productData.Date,productData.Equity,'equity','daily');

% ���� mPortfolio ��Ķ���
portfolio = mPortfolio(equity,posData(:,namesPosData),tradeData(:,namesTradeData));
%% ���� Instruments, positions �� trades

% ���� Instruments. create Instruments
% ����Ĳ���ΪInstrument������
portfolio.addInstruments('Multiplier')

% Ϊÿһ��Instrument����һ�� position ��
% create position data table for each symbol/Instrument
portfolio.cpositions()

% Ϊÿһ��Instrument����һ�� trade ��
% create trades data table for each symbol
portfolio.ctrades()

%% ��������������
% �����պ��������ʣ�Method��ѡ Simple �� Continuous ,Ĭ��Simple
% ��������˵Ӧ��ʹ�ö��������ʼ����׼�������
portfolio.calcDailyReturn('Simple') ;
portfolio.calcWeeklyReturn('Simple');
portfolio.calcMonthlyReturn('Simple');

%% �������ͳ����,������������̬�ֲ�

% �������ͳ����
portfolio.calcBasicStats('daily');
portfolio.calcBasicStats('weekly');
portfolio.calcBasicStats('monthly');

% �������������Ƿ������̬�ֲ�
testNormality(fts2mat(portfolio.returns.daily));




%% ���س�
portfolio.calcMaxDrawdown();
% ����س�ͳ������
portfolio.calcDrawdownStats();
% portfolio.drawdownStats

%% �������ձ���
% �趨�޷������ʣ����޷����ʲ����껯������
% ��ʹ��0, ���ڹ�ծ���ʻ�SHIBOR����
riskFreeRate = 0;
portfolio.calcSharpe(riskFreeRate,'daily'); %����������
portfolio.calcSharpe(riskFreeRate,'weekly'); %����������
portfolio.calcSharpe(riskFreeRate,'monthly'); %����������


%% Ratios

%����Sortino ratio

% Ŀ���������껯ֵ
MAR = 0;

portfolio.calcSortino( MAR, 'daily');
portfolio.calcSortino( MAR, 'weekly');
portfolio.calcSortino( MAR, 'monthly');

% ���� upside potential ratio
portfolio.calcUpsidePotential(MAR,'daily');
portfolio.calcUpsidePotential(MAR,'weekly');
portfolio.calcUpsidePotential(MAR,'monthly');

% ���� Sterling ratio
portfolio.calcSterling(MAR,'daily');
portfolio.calcSterling(MAR,'weekly');
portfolio.calcSterling(MAR,'monthly');

portfolio.calcSterling_maxDD(MAR,'daily');
portfolio.calcSterling_maxDD(MAR,'weekly');
portfolio.calcSterling_maxDD(MAR,'monthly');

% ���� Burke ratio
portfolio.calcBurke(MAR,'daily');
portfolio.calcBurke(MAR,'weekly');
portfolio.calcBurke(MAR,'monthly');

portfolio.calcBurke_maxDD(MAR,'daily');
portfolio.calcBurke_maxDD(MAR,'weekly');
portfolio.calcBurke_maxDD(MAR,'monthly');

% �µ�Ƶ��
portfolio.calcDownsideFreq(MAR, 'monthly');
portfolio.calcDownsideFreq(MAR, 'weekly');
portfolio.calcDownsideFreq(MAR, 'daily');
% ���������ʵ����ձ���
portfolio.calcSharpe_upside(MAR, 'monthly');
portfolio.calcSharpe_upside(MAR, 'weekly');
portfolio.calcSharpe_upside(MAR, 'daily');
% ���������ʵ����ձ���
portfolio.calcSharpe_downside(MAR, 'monthly' );
portfolio.calcSharpe_downside(MAR, 'weekly' );
portfolio.calcSharpe_downside(MAR, 'daily' );

%% ������ʳɽ���¼����ҵ������

% ����Ƿ��к�Լͬʱ���ж�ղ�
% �����ʾ����ͬһ��Լ�ϣ���ͷ�Ϳ�ͷ�Ľ���Ӧ���������м��㣬�Ի�ԭͬʱ��˫�ֲֵ߳����
portfolio.checkLongShortPos

% ���� trades �е�ÿ�����ڹ��ߵĽ������ݼ�����ʽ���ӯ��
portfolio.calcTxn_allInstruments

% ͳ�ƴӿ��ֵ�ƽ�ֵ�ÿ�ʽ���
portfolio.calcPerTradeStats
% ������Ʒ��ͳ��ӯ������
portfolio.calcAllSymbolsPL
% �������н���Ʒ�ֵ�����ӯ������
portfolio.aggregateSymbolPL

% ����ÿ�����׵ĺ�Լ��ÿ���λ��ӯ��
portfolio.calcPositionsPL

%% ��ʾ����
portfolio.dispAllSymbolPL

% �鿴������Լ�Ľ�������ͳ��
% ��ME501Ϊ��
portfolio.perTradeStats.ME505

%% ��ͼ

% �������ʷֲ�ͼ
portfolio.chartReturnDistribution('daily')

% �¶�������ֱ��ͼ
portfolio.chartMonthlyReturn

% chart position ������һ����Լ�ϵĲ�λ,ӯ������������

% ����Դ����ֻ�н���ۣ������������Խ���۳��ֵ�
% ��ME501Ϊ����
portfolio.chartPosition('ME505')
% ��ÿ���µ����ձ��ʵ�ͼ
portfolio.chartSharpe_eachmonth
% ��Ȩ������
portfolio.chartEquity
% ���س���ʱ������ͼ
portfolio.chartDrawdown


%% ����������ļ�
% ������н���Ʒ�ֵ�ӯ��
portfolio.outAllSymbolPL
% ����ͳ����
portfolio.outBasicStats

% ���ҵ�������ı���
portfolio.outRatios

% ������ϵ����еĽ���Ʒ�ֵ�ÿ�ʽ���ͳ������
portfolio.outPerTradeStats
% ����س�ͳ������
portfolio.outDrawdownStats

