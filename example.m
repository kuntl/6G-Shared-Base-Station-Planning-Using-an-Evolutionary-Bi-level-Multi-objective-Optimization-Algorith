tic
clear;
%�������ŵ�����
rb_down_num = 150;
rb_up_num = 80;

%�������cs_num����ѡվַ
cs_num = 80;
area = 4000;
CS = area*rand(cs_num,2);  %��ѡվַ
%��ÿ����ѡվַ�����վ�ķ���
r_tmp = rand(cs_num,1);
CS_COST = 30+50*r_tmp;

%��ѡվַ�س���·�ĵ�λ����
Bh_C_u = 0.00000001*r_tmp;

%ÿ����ѡվַ�����
CS_RENT = 1.5*CS_COST; %δ���ǹ����ۿ�

%ÿ����Ӫ���������tp_num�����Ե�
op_num = 3;
tp_num = 500;
TP =area*rand(op_num*tp_num,2);

%������˾�ں�ѡվַ�������ѡbs_num��λ�ý����վ
bs_num = cs_num;
%bs_num = randi([round(cs_num/2),cs_num]);
index_bs = randperm(cs_num,bs_num); %ѡ����վַ���  
BS = zeros(bs_num,2);     
%������˾�����վ�ĳɱ�
InC = sum(CS_COST(index_bs));
%������˾�����վ�ͻس���·�ĳɱ�
Bh_BS_u = Bh_C_u(index_bs);
BS_RENT = CS_RENT(index_bs);

% for i = 1:bs_num
%     BS(i,:) = CS(index_bs(i),:);
% end
BS = CS(index_bs,:);  %ѡ����վַ

%ÿ����Ӫ��ѡ�����û�վ����6G���磬����ÿ����Ӫ������20����վ
%op_rent_num = 20;
op_rent_num = randi([round(bs_num/2),bs_num],op_num,1);
op_rent_bs = zeros(sum(op_rent_num),2);
%op_rent_bs = zeros(op_rent_num*op_num,2);  %��Ӫ�����õ�վַ
op_frequency = [230,250,270,290,310,330];

%ͳ�ƻ�վ����Ӫ�̵�����
bs_rent_count = zeros(bs_num,1);  %��վ����Ӫ�̵�����

%ÿ����Ӫ�����õĻ�վ���
%index_ops_bs = zeros(op_num,op_rent_num);
index_ops_bs = zeros(sum(op_rent_num),1);
for i = 1:1:op_num
    if i >=2
        tmp_num = op_rent_num(1:(i-1),:);
        row_bs = (sum(tmp_num(:))+1):(sum(tmp_num(:))+op_rent_num(i));%��Ӫ��i�Ļ�վ����
    else
        row_bs = 1:op_rent_num(i);%��Ӫ��i�Ļ�վ����
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

%ͳ�ƻ�վ����Ӫ�̵�����
%bs_rent_count = zeros(bs_num,1);  %��վ����Ӫ�̵�����

%iob_tmp = unique(index_ops_bs);
%bs_rent_count(iob_tmp) = histc(index_ops_bs(:),iob_tmp);

%ÿ����վ�Ĺ����ۿ�
D = 1+0.2*bs_rent_count;
BS_RENT = BS_RENT./D;

%��Ӫ����Ҫ֧�������
% ops_rent = zeros(op_num,1);
% for i = 1:1:op_num
%     index_op_bs = index_ops_bs(i,:);
%     ops_rent(i) = sum(BS_RENT(index_op_bs));    
% end
ops_rent = 0;
for i = 1:1:op_num
    if i >=2
        tmp_num = op_rent_num(1:(i-1),:);
        row_bs = (sum(tmp_num(:))+1):(sum(tmp_num(:))+op_rent_num(i));%��Ӫ��i�Ļ�վ����
    else
        row_bs = 1:op_rent_num(i);%��Ӫ��i�Ļ�վ����
    end
    ops_rent = ops_rent + sum(BS_RENT(index_ops_bs(row_bs)));
end


%ȷ���������ò���
height_type = 4;
total_power_type = 10;
op_rent_bs_ac(:,1) = randi(height_type,sum(op_rent_num),1);        %���߸߶�
op_rent_bs_ac(:,2) = randi(total_power_type,sum(op_rent_num),1);   %���߹���

op_ac_cost = zeros(op_num,1);
P_tot = 0;
for i = 1:1:op_num
    if i >=2
        tmp_num = op_rent_num(1:(i-1),:);
        row_bs = (sum(tmp_num(:))+1):(sum(tmp_num(:))+op_rent_num(i));%��Ӫ��i�Ļ�վ����
    else
        row_bs = 1:op_rent_num(i);%��Ӫ��i�Ļ�վ����
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


%�û����������������
rate_up_min = 1e4;
rate_down_min = 1e5;

%��װ����
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
%������Ե������������
[tps_down_rate,tps_up_rate,tp_bs_allo] = calRateofTps_v3(TP,op_rent_bs,op_rent_bs_ac,parameter);

%��װ����
parameter.index_ops_bs = index_ops_bs; %ÿ����Ӫ�����û�վ�ı��
parameter.Bh_BS_u = Bh_BS_u; %ÿ����վ��λ�س���·����ɱ�

BhC = calBhCofBSs(tps_down_rate,tps_up_rate,tp_bs_allo,parameter);

cov = zeros(tp_num*op_num,1);
[cov_size,~] = size(cov);
%���㸲���ʺ�ϵͳ����
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

%������˾�����վ�ĳɱ�
InC
%������˾����س���·�ĳɱ�
BhC_sum = sum(BhC)
%������˾����Դ����
P_tot
%���������ڵ�ų���¶
EMF
%���縲����
fcov = sum(cov)/(tp_num*op_num) 
%��������
fcap
%��Ӫ�̵����ɱ�
ops_rent_sum = sum(ops_rent)
%��Ӫ�̵����߽���ɱ�
AC_sum = sum(op_ac_cost)

toc
%������վ�滮ʾ��ͼ
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




