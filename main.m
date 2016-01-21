%% 说明
% 1. MATLAB 版本问题(代码基于MATLAB R2014b版)：
% 1.1 MATLAB引入了table数据类型用于替换dataset类型, 在以下代码中使用了table
%     如果由于MATLAB版本问题无法调用table，可使用功能相同的dataset
% 1.2 fts (financial time series)类型需要 Financial toolbox 的支持
%
% 2. mPortfolio 类
%   根据面向对象程序设计的思想自定义并使用了 mPortfolio 类（前缀m用于和Financial toolbox的portfolio类区别)
%   尽量将对投资组合所进行的计算和统计的函数封装进了 mPortfolio 类中，可在mPortfolio.m文件中查看修改
%   使用类封装可以解决变量名冲突，使对象结构清晰，但是在注重效率的应用中，时间开销是需要考虑的问题
%   
% 3. 运行代码： 请按section顺序依次运行代码导入并创建portfolio对象，业绩分析相关代码可重复运行
%% 导入数据集
% 导入数据集
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

%% 设置导入的数据表列名称
% 以下命令可查看导入后的数据集的列名称
% these commands can access the variable names of the tables; 
% VariableNames need to be set manually.
% (dataset)
% posData.Properties.VarNames
% productData.Properties.VarNames
% tradeData.Properties.VarNames
% (table)
% 
% 
% 该命令可用于查询数据集各变量的数据类型
% (dataset)
% datasetfun(@class,posData,'UniformOutput',false)
% (table)
% posData.Properties.VariableDescriptions

% modify the names manually
% 修改各列的名称
namesPosData = {'Date','Symbol',...
    'LongPos','ShortPos','PrevStlPrice','StlPrice','Multiplier'};
posData.Properties.VariableNames([1,2,3,5,7,8,9]) = namesPosData; 
namesProductData = {'Date','Equity'};
productData.Properties.VariableNames([1,2]) = namesProductData;
namesTradeData = {'Date',...
    'Symbol','LongQty','LongPrice','ShortQty','ShortPrice','PrevStlPrice',...
    'StlPrice','Multiplier','TxnFee'};
tradeData.Properties.VariableNames([1,2,3,4,5,6,7,8,9,10]) = namesTradeData;


%% clean the Nan row in tables. 针对导入的数据表进行空行清理
% 清理数据表中Date为空的行
productData = cleanNanRow(productData,'Date');
posData = cleanNanRow(posData,'Date');
tradeData = cleanNanRow(tradeData,'Date');

% 清理 tradeData 中没有操作记录的行. cellstring中空值等于 ''
tradeData = cleanNanRow(tradeData,'Symbol');

%% 针对导入的数据表进行数据类型转换
%将datestring转换为 datenum类型
posData.Date = datenum(posData.Date);
productData.Date = datenum(num2str(productData.Date),'yyyymmdd');
tradeData.Date = datenum( tradeData.Date);

% 将Symbol由cell array of string转换为categorical类型
posData.Symbol = categorical(posData.Symbol);
tradeData.Symbol = categorical(tradeData.Symbol);
%% 创建portfolio对象

% 创建权益的时间序列. create financial time series object of equity
equity = fints(productData.Date,productData.Equity,'equity','daily');

% 创建 mPortfolio 类的对象
portfolio = mPortfolio(equity,posData(:,namesPosData),tradeData(:,namesTradeData));
%% 创建 Instruments, positions 和 trades

% 创建 Instruments. create Instruments
% 输入的参数为Instrument的属性
portfolio.addInstruments('Multiplier')

% 为每一个Instrument创建一个 position 表
% create position data table for each symbol/Instrument
portfolio.cpositions()

% 为每一个Instrument创建一个 trade 表
% create trades data table for each symbol
portfolio.ctrades()

%% 计算收益率序列
% 计算日和月收益率，Method可选 Simple 或 Continuous ,默认Simple
% 理论上来说应该使用对数收益率计算标准差即波动率
portfolio.calcDailyReturn('Simple') ;
portfolio.calcWeeklyReturn('Simple');
portfolio.calcMonthlyReturn('Simple');

%% 计算基本统计量,检验收益率正态分布

% 计算基本统计量
portfolio.calcBasicStats('daily');
portfolio.calcBasicStats('weekly');
portfolio.calcBasicStats('monthly');

% 检验日收益率是否符合正态分布
testNormality(fts2mat(portfolio.returns.daily));




%% 最大回撤
portfolio.calcMaxDrawdown();
% 计算回撤统计数据
portfolio.calcDrawdownStats();
% portfolio.drawdownStats

%% 计算夏普比率
% 设定无风险利率，即无风险资产的年化收益率
% 可使用0, 短期国债利率或SHIBOR利率
riskFreeRate = 0;
portfolio.calcSharpe(riskFreeRate,'daily'); %基于日收益
portfolio.calcSharpe(riskFreeRate,'weekly'); %基于周收益
portfolio.calcSharpe(riskFreeRate,'monthly'); %基于月收益


%% Ratios

%计算Sortino ratio

% 目标收益率年化值
MAR = 0;

portfolio.calcSortino( MAR, 'daily');
portfolio.calcSortino( MAR, 'weekly');
portfolio.calcSortino( MAR, 'monthly');

% 计算 upside potential ratio
portfolio.calcUpsidePotential(MAR,'daily');
portfolio.calcUpsidePotential(MAR,'weekly');
portfolio.calcUpsidePotential(MAR,'monthly');

% 计算 Sterling ratio
portfolio.calcSterling(MAR,'daily');
portfolio.calcSterling(MAR,'weekly');
portfolio.calcSterling(MAR,'monthly');

portfolio.calcSterling_maxDD(MAR,'daily');
portfolio.calcSterling_maxDD(MAR,'weekly');
portfolio.calcSterling_maxDD(MAR,'monthly');

% 计算 Burke ratio
portfolio.calcBurke(MAR,'daily');
portfolio.calcBurke(MAR,'weekly');
portfolio.calcBurke(MAR,'monthly');

portfolio.calcBurke_maxDD(MAR,'daily');
portfolio.calcBurke_maxDD(MAR,'weekly');
portfolio.calcBurke_maxDD(MAR,'monthly');

% 下跌频率
portfolio.calcDownsideFreq(MAR, 'monthly');
portfolio.calcDownsideFreq(MAR, 'weekly');
portfolio.calcDownsideFreq(MAR, 'daily');
% 超额收益率的夏普比率
portfolio.calcSharpe_upside(MAR, 'monthly');
portfolio.calcSharpe_upside(MAR, 'weekly');
portfolio.calcSharpe_upside(MAR, 'daily');
% 下行收益率的夏普比率
portfolio.calcSharpe_downside(MAR, 'monthly' );
portfolio.calcSharpe_downside(MAR, 'weekly' );
portfolio.calcSharpe_downside(MAR, 'daily' );

%% 根据逐笔成交记录进行业绩分析

% 检查是否有合约同时持有多空仓
% 检查显示，在同一合约上，多头和空头的交易应独立开进行计算，以还原同时在双边持仓的情况
portfolio.checkLongShortPos

% 根据 trades 中的每个金融工具的交易数据计算逐笔交易盈亏
portfolio.calcTxn_allInstruments

% 统计从开仓到平仓的每笔交易
portfolio.calcPerTradeStats
% 按交易品种统计盈亏数据
portfolio.calcAllSymbolsPL
% 整合所有交易品种的整体盈亏数据
portfolio.aggregateSymbolPL

% 计算每个交易的合约的每天仓位的盈亏
portfolio.calcPositionsPL

%% 显示数据
portfolio.dispAllSymbolPL

% 查看单个合约的交易数据统计
% 以ME501为例
portfolio.perTradeStats.ME505

%% 绘图

% 日收益率分布图
portfolio.chartReturnDistribution('daily')

% 月度收益率直方图
portfolio.chartMonthlyReturn

% chart position 画出在一个合约上的仓位,盈亏和行情曲线

% 由于源数据只有结算价，行情曲线是以结算价呈现的
% 以ME501为例：
portfolio.chartPosition('ME505')
% 画每个月的夏普比率的图
portfolio.chartSharpe_eachmonth
% 画权益曲线
portfolio.chartEquity
% 画回撤的时间序列图
portfolio.chartDrawdown


%% 数据输出到文件
% 输出所有交易品种的盈亏
portfolio.outAllSymbolPL
% 基本统计量
portfolio.outBasicStats

% 输出业绩分析的比率
portfolio.outRatios

% 输出整合的所有的交易品种的每笔交易统计数据
portfolio.outPerTradeStats
% 输出回撤统计数据
portfolio.outDrawdownStats

