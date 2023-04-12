function BhC = calBhCofBSs(tps_down_rate,tps_up_rate,tp_bs_allo,parameter)

op_num = parameter.op_num;  
tp_num = parameter.tp_num;   
bs_num= parameter.bs_num;    
index_ops_bs = parameter.index_ops_bs ; 
Bh_BS_u = parameter.Bh_BS_u ; 
op_rent_num = parameter.op_rent_num;
BhC = zeros(bs_num,1); 

for i = 1:1:op_num
    if op_rent_num(i)>0
    row_tp = ((i-1)*tp_num+1):((i-1)*tp_num+tp_num);
    if i >=2
        tmp_num = op_rent_num(1:(i-1),:);
        row_bs = (sum(tmp_num(:))+1):(sum(tmp_num(:))+op_rent_num(i));
    else
        row_bs = 1:op_rent_num(i);
    end
    index_op_bs = index_ops_bs(row_bs,:);
    tp_bs_allo_op = tp_bs_allo(row_tp);
    tps_down_rate_op = tps_down_rate(row_tp);
    tps_up_rate_op = tps_up_rate(row_tp);
    for u = 1:1:tp_num
        bs_index = index_op_bs(tp_bs_allo_op(u));
        BhC(bs_index) = BhC(bs_index)+Bh_BS_u(bs_index)*(tps_down_rate_op(u)+tps_up_rate_op(u));
    end    
    end
end

%%
% Row_tp = reshape(1:op_num*tp_num,tp_num,1,op_num);
% tp_bs_allo_ops = tp_bs_allo(Row_tp);
% tps_down_rate_ops = tps_down_rate(Row_tp);
% tps_up_rate_ops = tps_up_rate(Row_tp);
% 
% Index_ops_bs = index_ops_bs(:);
% I = triu(ones(op_num,op_num),0);
% I = [zeros(op_num,1),I(:,1:op_num-1)];
% bs_index = Index_ops_bs(tp_bs_allo_ops+reshape(sum(op_rent_num.*I),1,1,op_num));
% 
% 
% if gpuDeviceCount>0 
%     Bh_BS_u = gpuArray(Bh_BS_u);
%     bs_index = gpuArray(bs_index);
%     tps_down_rate_ops = gpuArray(tps_down_rate_ops);
%     tps_up_rate_ops = gpuArray(tps_up_rate_ops);
% end
% %%
% BhCs = Bh_BS_u(bs_index).*(tps_down_rate_ops+tps_up_rate_ops);
% site = bs_index(:)==[1:bs_num];
% BhC = (sum(repmat(BhCs(:),1,bs_num).*site))';
% 
% if gpuDeviceCount>0 
%     BhC = gather(BhC);
% end
end