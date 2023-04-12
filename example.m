tic
clear;
%上下行信道数量
rb_down_num = 150;
rb_up_num = 80;

%随机生成cs_num个候选站址
cs_num = 80;
area = 4000;
CS = area*rand(cs_num,2);  %候选站址
%在每个候选站址建设基站的费用
r_tmp = rand(cs_num,1);
CS_COST = 30+50*r_tmp;

%候选站址回程链路的单位费用
Bh_C_u = 0.00000001*r_tmp;

%每个候选站址的租金
CS_RENT = 1.5*CS_COST; %未考虑共享折扣

%每个运营商随机生成tp_num个测试点
op_num = 3;
tp_num = 500;
TP =area*rand(op_num*tp_num,2);

%铁塔公司在候选站址中随机挑选bs_num个位置建设基站
bs_num = cs_num;
%bs_num = randi([round(cs_num/2),cs_num]);
index_bs = randperm(cs_num,bs_num); %选出的站址编号  
BS = zeros(bs_num,2);     
%铁塔公司建设基站的成本
InC = sum(CS_COST(index_bs));
%铁塔公司建设基站和回程链路的成本
Bh_BS_u = Bh_C_u(index_bs);
BS_RENT = CS_RENT(index_bs);

% for i = 1:bs_num
%     BS(i,:) = CS(index_bs(i),:);
% end
BS = CS(index_bs,:);  %选出的站址

%每个运营商选择租用基站建设6G网络，假设每个运营商租用20个基站
%op_rent_num = 20;
op_rent_num = randi([round(bs_num/2),bs_num],op_num,1);
op_rent_bs = zeros(sum(op_rent_num),2);
%op_rent_bs = zeros(op_rent_num*op_num,2);  %运营商租用的站址
op_frequency = [230,250,270,290,310,330];

%统计基站上运营商的数量
bs_rent_count = zeros(bs_num,1);  %基站上运营商的数量

%每个运营商租用的基站编号
%index_ops_bs = zeros(op_num,op_rent_num);
index_ops_bs = zeros(sum(op_rent_num),1);
for i = 1:1:op_num
    if i >=2
        tmp_num = op_rent_num(1:(i-1),:);
        row_bs = (sum(tmp_num(:))+1):(sum(tmp_num(:))+op_rent_num(i));%运营商i的基站行数
    else
        row_bs = 1:op_rent_num(i);%运营商i的基站行数
    end
    index_ops_bs(row_bs) = randperm(bs_num,op_rent_num(i));
    bs_rent_count(index_ops_bs(row_bs)) = bs_rent_count(index_ops_bs(row_bs))+1;
%     k = 1;
%     for j = ((i-1)*op_rent_num+1):((i-1)*op_rent_num+op_rent_num)
%         op_rent_bs(j,:) = BS(index_op_bs(k),:);
%         k = k+1;
%     end
    op_rent_bs(row_bs,:) = BS(index_ops_bs(row_bs),:);
end
% Idx_ops_bs = index_ops_bs';
% op_rent_bs = BS(Idx_ops_bs(:),:);

%统计基站上运营商的数量
%bs_rent_count = zeros(bs_num,1);  %基站上运营商的数量

%iob_tmp = unique(index_ops_bs);
%bs_rent_count(iob_tmp) = histc(index_ops_bs(:),iob_tmp);

%每个基站的共享折扣
D = 1+0.2*bs_rent_count;
BS_RENT = BS_RENT./D;

%运营商需要支付的租金
% ops_rent = zeros(op_num,1);
% for i = 1:1:op_num
%     index_op_bs = index_ops_bs(i,:);
%     ops_rent(i) = sum(BS_RENT(index_op_bs));    
% end
ops_rent = 0;
for i = 1:1:op_num
    if i >=2
        tmp_num = op_rent_num(1:(i-1),:);
        row_bs = (sum(tmp_num(:))+1):(sum(tmp_num(:))+op_rent_num(i));%运营商i的基站行数
    else
        row_bs = 1:op_rent_num(i);%运营商i的基站行数
    end
    ops_rent = ops_rent + sum(BS_RENT(index_ops_bs(row_bs)));
end


%确定天线配置参数
height_type = 4;
total_power_type = 10;
op_rent_bs_ac(:,1) = randi(height_type,sum(op_rent_num),1);        %天线高度
op_rent_bs_ac(:,2) = randi(total_power_type,sum(op_rent_num),1);   %天线功率

op_ac_cost = zeros(op_num,1);
P_tot = 0;
for i = 1:1:op_num
    if i >=2
        tmp_num = op_rent_num(1:(i-1),:);
        row_bs = (sum(tmp_num(:))+1):(sum(tmp_num(:))+op_rent_num(i));%运营商i的基站行数
    else
        row_bs = 1:op_rent_num(i);%运营商i的基站行数
    end

    op_rent_bs_ac_i = op_rent_bs_ac(row_bs,:);
%     for j = 1:1:op_rent_num
%         [total_power,~,cost] = ac_decode(op_rent_bs_ac_i(j,:));
%         op_ac_cost(i) = op_ac_cost(i) + cost; 
%         P_tot = P_tot + total_power;
%     end
    [total_power,~,cost] = ac_decode(op_rent_bs_ac_i);
    op_ac_cost(i) = sum(cost);
    P_tot = P_tot+sum(total_power);
end


%用户的最低上下行速率
rate_up_min = 1e4;
rate_down_min = 1e5;

%封装参数
parameter.rb_down_num = rb_down_num;
parameter.rb_up_num = rb_up_num;
parameter.cs_num = cs_num;
parameter.area = area;
parameter.op_num = op_num;
parameter.tp_num = tp_num;
parameter.bs_num = bs_num;

parameter.op_rent_num = op_rent_num;

parameter.op_frequency = op_frequency;
parameter.rate_up_min = rate_up_min;
parameter.rate_down_min = rate_down_min;
parameter.area = area;

EMF = calEMFofArea(op_rent_bs,op_rent_bs_ac,parameter);
%计算测试点的上下行速率
[tps_down_rate,tps_up_rate,tp_bs_allo] = calRateofTps_v3(TP,op_rent_bs,op_rent_bs_ac,parameter);

%封装参数
parameter.index_ops_bs = index_ops_bs; %每个运营商租用基站的编号
parameter.Bh_BS_u = Bh_BS_u; %每个基站单位回程链路建设成本

BhC = calBhCofBSs(tps_down_rate,tps_up_rate,tp_bs_allo,parameter);

cov = zeros(tp_num*op_num,1);
[cov_size,~] = size(cov);
%计算覆盖率和系统容量
%sum(tps_down_rate)
%sum(tps_up_rate)
fcap = sum(tps_down_rate)+sum(tps_up_rate);
% for u = 1:1:cov_size
%     if tps_down_rate(u) > rate_down_min && tps_up_rate(u) > rate_up_min
%         cov(u) = 1;
%     end
%     
% end
cov = tps_down_rate > rate_down_min & tps_up_rate > rate_up_min;

%铁塔公司建设基站的成本
InC
%铁塔公司建设回程链路的成本
BhC_sum = sum(BhC)
%铁塔公司的能源消耗
P_tot
%服务区域内电磁场暴露
EMF
%网络覆盖率
fcov = sum(cov)/(tp_num*op_num) 
%网络容量
fcap
%运营商的租金成本
ops_rent_sum = sum(ops_rent)
%运营商的天线建设成本
AC_sum = sum(op_ac_cost)

toc
%画出基站规划示意图
% figure(1);
% plot(CS(:,1)',CS(:,2),'ok','MarkerSize',5);
% hold on;
% plot(BS(:,1)',BS(:,2),'or','MarkerSize',5);
% hold on;
% for i = 1:1:op_num
%     figure(i+1);
%     plot(BS(:,1)',BS(:,2),'ok','MarkerSize',5);
%     hold on;
%     row_tp = ((i-1)*tp_num+1):((i-1)*tp_num+tp_num);
%     TP_op = TP(row_tp,:);
%     plot(TP_op(:,1)',TP_op(:,2),'xr','MarkerSize',8);
%     hold on;
%     row_bs = ((i-1)*op_rent_num+1):((i-1)*op_rent_num+op_rent_num);
%     op_rent_bs_i = op_rent_bs(row_bs,:);
%     plot(op_rent_bs_i(:,1)',op_rent_bs_i(:,2),'or','MarkerSize',5);
%     hold on;
%     for u = 1:1:tp_num
%         if cov(u+(i-1)*tp_num) == 0
%             plot(TP(u,1),TP(u,2),'xk','MarkerSize',8);
%             hold on;
%         end
%     end
% end




