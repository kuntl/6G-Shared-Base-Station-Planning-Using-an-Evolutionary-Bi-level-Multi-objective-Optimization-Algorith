function [tps_down_rate,tps_up_rate,tps_bs_allo]=calRateofTps_v2(TP,op_rent_bs,op_rent_bs_ac,parameter)

rb_down_num = parameter.rb_down_num;
rb_up_num = parameter.rb_up_num;
op_num = parameter.op_num;
tp_num = parameter.tp_num;
op_rent_num = parameter.op_rent_num;
op_frequency = parameter.op_frequency;
rate_up_min = parameter.rate_up_min;
rate_down_min = parameter.rate_down_min;
tps_bs_allo = zeros(tp_num*op_num,1);
tps_down_rate = zeros(tp_num*op_num,1);
tps_up_rate = zeros(tp_num*op_num,1);

for i = 1:1:op_num
    op_rent_num_i = op_rent_num(i);
    row_tp = ((i-1)*tp_num+1):((i-1)*tp_num+tp_num);
    if op_rent_num_i >0
        %row_bs = ((i-1)*op_rent_num+1):((i-1)*op_rent_num+op_rent_num);
        if i >=2
            tmp_num = op_rent_num(1:(i-1),:);
            row_bs = (sum(tmp_num(:))+1):(sum(tmp_num(:))+op_rent_num_i);
        else
            row_bs = 1:op_rent_num_i;
        end
        
        TP_op = TP(row_tp,:);
        op_rent_bs_i = op_rent_bs(row_bs,:);
        op_frequency_i = op_frequency(i);
        op_rent_bs_ac_i = op_rent_bs_ac(row_bs,:);
        [total_power,height,~] = ac_decode(op_rent_bs_ac_i);
        
        %% 
        %     if gpuDeviceCount>0
        %
        %         height = gpuArray(height);
        %         total_power = gpuArray(total_power);
        %         tp_op = gpuArray(tp_op);
        %
        %         op_rent_bs_i = gpuArray(op_rent_bs_i);
        %
        %     %%
        %
        %     Tp_bs_distance = sqrt(sum([repmat(height,1,1,tp_num)-1.5,(tp_op-op_rent_bs_i)].^2,2));
        %     switch op_frequency_i
        %         case 230
        %             op_gain_down = (1/10).^((81.10+21.61*log10(Tp_bs_distance)+(0+5.29.*gpuArray.randn(op_rent_num_i,rb_down_num,tp_num))-20)/10);
        %         case 250
        %             op_gain_down = (1/10).^((83.21+20.09*log10(Tp_bs_distance)+(0+8.41.*gpuArray.randn(op_rent_num_i,rb_down_num,tp_num))-20)/10);
        %         case 270
        %             op_gain_down = (1/10).^((83.24+18.21*log10(Tp_bs_distance)+(0+9.61.*gpuArray.randn(op_rent_num_i,rb_down_num,tp_num))-20)/10);
        %         case 290
        %             op_gain_down = (1/10).^((85.35+18.31*log10(Tp_bs_distance)+(0+10.24.*gpuArray.randn(op_rent_num_i,rb_down_num,tp_num))-20)/10);
        %         case 310
        %             op_gain_down = (1/10).^((80.56+17.51*log10(Tp_bs_distance)+(0+7.84.*gpuArray.randn(op_rent_num_i,rb_down_num,tp_num))-20)/10);
        %         case 330
        %             op_gain_down = (1/10).^((85.02+17.45*log10(Tp_bs_distance)+(0+9.*gpuArray.randn(op_rent_num_i,rb_down_num,tp_num))-20)/10);
        %     end
        %     switch op_frequency_i
        %         case 230
        %             op_gain_up = (1/10).^((81.10+21.61*log10(Tp_bs_distance)+(0+5.29.*gpuArray.randn(op_rent_num_i,rb_up_num,tp_num))-20)/10);
        %         case 250
        %             op_gain_up = (1/10).^((83.21+20.09*log10(Tp_bs_distance)+(0+8.41.*gpuArray.randn(op_rent_num_i,rb_up_num,tp_num))-20)/10);
        %         case 270
        %             op_gain_up = (1/10).^((83.24+18.21*log10(Tp_bs_distance)+(0+9.61.*gpuArray.randn(op_rent_num_i,rb_up_num,tp_num))-20)/10);
        %         case 290
        %             op_gain_up = (1/10).^((85.35+18.31*log10(Tp_bs_distance)+(0+10.24.*gpuArray.randn(op_rent_num_i,rb_up_num,tp_num))-20)/10);
        %         case 310
        %             op_gain_up = (1/10).^((80.56+17.51*log10(Tp_bs_distance)+(0+7.84.*gpuArray.randn(op_rent_num_i,rb_up_num,tp_num))-20)/10);
        %         case 330
        %             op_gain_up = (1/10).^((85.02+17.45*log10(Tp_bs_distance)+(0+9.*gpuArray.randn(op_rent_num_i,rb_up_num,tp_num))-20)/10);
        %     end
        %     tp_bs_v_gain = (reshape(total_power.*mean(op_gain_down,2),op_rent_num_i,tp_num))';
        %
        %       [~,tp_bs_allo_i] = max(tp_bs_v_gain,[],2); 
        %
        %     
        %
        %         tp_bs_allo_i = gather(tp_bs_allo_i);
        %         op_gain_down = gather(op_gain_down);
        %         op_gain_up = gather(op_gain_up);
        %         total_power = gather(total_power);
        %
        %         clear height op_rent_bs_i tp_bs_v_gain Tp_bs_distance
        %     else
        
        % Tp_bs_distance = sqrt(sum([repmat(height,1,1,tp_num)-1.5,(tp_op-op_rent_bs_i)].^2,2));
        %tp_op = reshape(TP_op,1,2,tp_num); 
        %op_rent_bs_ = reshape(op_rent_bs(row_bs,:),1,2,op_rent_num_i);
        %tmpp =  TP_op - op_rent_bs_; 
        %tmp = tp_op-op_rent_bs_i;
        %dis = repmat(height,1,1,tp_num);
        %diss = reshape(height,1,1,op_rent_num_i)+zeros(300,1);
        %tmp_2 = [repmat(height,1,1,tp_num)-1.5,reshape(TP_op,1,2,tp_num)-op_rent_bs_i].^2;
        %tmp_22 = [diss-1.5, tmpp].^2;
        %Tp_bs_distance_ = sqrt(sum(tmp_22,2));
        
        
        Tp_bs_distance = sqrt(sum([repmat(height,1,1,tp_num)-1.5,reshape(TP_op,1,2,tp_num)-op_rent_bs_i].^2,2));
        
        [value, index_]=sort(Tp_bs_distance,1,'descend');
        
%         if op_rent_num_i > 5
%             s_num = 5;
%         else
%             s_num = op_rent_num_i;
%         end
        s_num = min(5,op_rent_num_i);
        
        Tp_bs_distance_index = index_(1:s_num,:,:);
        
        Tp_bs_distance_ = value(1:s_num,:,:);
       
      if gpuDeviceCount>0 
          Tp_bs_distance_ = gpuArray(single(Tp_bs_distance_));
           
             switch op_frequency_i
                case 230
                    op_gain_down = (1/10).^((81.10+21.61*log10(Tp_bs_distance_)+(0+5.29.*gpuArray.randn(s_num,rb_down_num,tp_num))-20)/10);
                case 250
                    op_gain_down = (1/10).^((83.21+20.09*log10(Tp_bs_distance_)+(0+8.41.*gpuArray.randn(s_num,rb_down_num,tp_num))-20)/10);
                case 270
                    op_gain_down = (1/10).^((83.24+18.21*log10(Tp_bs_distance_)+(0+9.61.*gpuArray.randn(s_num,rb_down_num,tp_num))-20)/10);
                case 290
                    op_gain_down = (1/10).^((85.35+18.31*log10(Tp_bs_distance_)+(0+10.24.*gpuArray.randn(s_num,rb_down_num,tp_num))-20)/10);
                case 310
                    op_gain_down = (1/10).^((80.56+17.51*log10(Tp_bs_distance_)+(0+7.84.*gpuArray.randn(s_num,rb_down_num,tp_num))-20)/10);
                case 330
                    op_gain_down = (1/10).^((85.02+17.45*log10(Tp_bs_distance_)+(0+9.*gpuArray.randn(s_num,rb_down_num,tp_num))-20)/10);
             end
            switch op_frequency_i
                case 230
                    op_gain_up = (1/10).^((81.10+21.61*log10(Tp_bs_distance_)+(0+5.29.*gpuArray.randn(s_num,rb_up_num,tp_num))-20)/10);
                case 250
                    op_gain_up = (1/10).^((83.21+20.09*log10(Tp_bs_distance_)+(0+8.41.*gpuArray.randn(s_num,rb_up_num,tp_num))-20)/10);
                case 270
                    op_gain_up = (1/10).^((83.24+18.21*log10(Tp_bs_distance_)+(0+9.61.*gpuArray.randn(s_num,rb_up_num,tp_num))-20)/10);
                case 290
                    op_gain_up = (1/10).^((85.35+18.31*log10(Tp_bs_distance_)+(0+10.24.*gpuArray.randn(s_num,rb_up_num,tp_num))-20)/10);
                case 310
                    op_gain_up = (1/10).^((80.56+17.51*log10(Tp_bs_distance_)+(0+7.84.*gpuArray.randn(s_num,rb_up_num,tp_num))-20)/10);
                case 330
                    op_gain_up = (1/10).^((85.02+17.45*log10(Tp_bs_distance_)+(0+9.*gpuArray.randn(s_num,rb_up_num,tp_num))-20)/10);
            end
          
          op_gain_down = gather(op_gain_down);
          op_gain_up = gather(op_gain_up);
      else
            switch op_frequency_i
                case 230
                    op_gain_down = (1/10).^((81.10+21.61*log10(Tp_bs_distance_)+(0+5.29.*randn(s_num,rb_down_num,tp_num))-20)/10);
                case 250
                    op_gain_down = (1/10).^((83.21+20.09*log10(Tp_bs_distance_)+(0+8.41.*randn(s_num,rb_down_num,tp_num))-20)/10);
                case 270
                    op_gain_down = (1/10).^((83.24+18.21*log10(Tp_bs_distance_)+(0+9.61.*gpuArray.randn(s_num,rb_down_num,tp_num))-20)/10);
                case 290
                    op_gain_down = (1/10).^((85.35+18.31*log10(Tp_bs_distance_)+(0+10.24.*randn(s_num,rb_down_num,tp_num))-20)/10);
                case 310
                    op_gain_down = (1/10).^((80.56+17.51*log10(Tp_bs_distance_)+(0+7.84.*randn(s_num,rb_down_num,tp_num))-20)/10);
                case 330
                    op_gain_down = (1/10).^((85.02+17.45*log10(Tp_bs_distance_)+(0+9.*randn(s_num,rb_down_num,tp_num))-20)/10);
             end
            switch op_frequency_i
                case 230
                    op_gain_up = (1/10).^((81.10+21.61*log10(Tp_bs_distance_)+(0+5.29.*randn(s_num,rb_up_num,tp_num))-20)/10);
                case 250
                    op_gain_up = (1/10).^((83.21+20.09*log10(Tp_bs_distance_)+(0+8.41.*randn(s_num,rb_up_num,tp_num))-20)/10);
                case 270
                    op_gain_up = (1/10).^((83.24+18.21*log10(Tp_bs_distance_)+(0+9.61.*gpuArray.randn(s_num,rb_up_num,tp_num))-20)/10);
                case 290
                    op_gain_up = (1/10).^((85.35+18.31*log10(Tp_bs_distance_)+(0+10.24.*gpuArray.randn(s_num,rb_up_num,tp_num))-20)/10);
                case 310
                    op_gain_up = (1/10).^((80.56+17.51*log10(Tp_bs_distance_)+(0+7.84.*randn(s_num,rb_up_num,tp_num))-20)/10);
                case 330
                    op_gain_up = (1/10).^((85.02+17.45*log10(Tp_bs_distance_)+(0+9.*randn(s_num,rb_up_num,tp_num))-20)/10);
            end
      end
      
        %         op_gain_down = (1/10).^((28+22*log10(Tp_bs_distance)+20*log10(op_frequency_i)+normrnd(0,8,[op_rent_num_i,rb_down_num,tp_num])-20)/10);
        %         op_gain_up = (1/10).^((28+22*log10(Tp_bs_distance)+20*log10(op_frequency_i)+normrnd(0,8,[op_rent_num_i,rb_up_num,tp_num])-20)/10);
        %     
        %     op_gain_down = (1/10).^((28+22*log10(Tp_bs_distance)+20*log10*(op_frequency_i)+normrnd(0,8,[op_rent_num_i,rb_down_num,tp_num])-20)/10);
        %
        %    
        %     op_gain_up = (1/10).^((28+22*log10(Tp_bs_distance)+20*log10*(op_frequency_i)+normrnd(0,8,[op_rent_num_i,rb_up_num,tp_num])-20)/10);
        
        %     switch op_frequency_i
        %         case 230
        %             op_gain_down = (1/10).^((81.10+21.61*log10(Tp_bs_distance)+normrnd(0,5.29,[op_rent_num_i,rb_down_num,tp_num])-20)/10);
        %         case 250
        %             op_gain_down = (1/10).^((83.21+20.09*log10(Tp_bs_distance)+normrnd(0,8.41,[op_rent_num_i,rb_down_num,tp_num])-20)/10);
        %         case 270
        %             op_gain_down = (1/10).^((83.24+18.21*log10(Tp_bs_distance)+normrnd(0,9.61,[op_rent_num_i,rb_down_num,tp_num])-20)/10);
        %         case 290
        %             op_gain_down = (1/10).^((85.35+18.31*log10(Tp_bs_distance)+normrnd(0,10.24,[op_rent_num_i,rb_down_num,tp_num])-20)/10);
        %         case 310
        %             op_gain_down = (1/10).^((80.56+17.51*log10(Tp_bs_distance)+normrnd(0,7.84,[op_rent_num_i,rb_down_num,tp_num])-20)/10);
        %         case 330
        %             op_gain_down = (1/10).^((85.02+17.45*log10(Tp_bs_distance)+normrnd(0,9,[op_rent_num_i,rb_down_num,tp_num])-20)/10);
        %     end
        %
        %   
        %     switch op_frequency_i
        %         case 230
        %             op_gain_up = (1/10).^((81.10+21.61*log10(Tp_bs_distance)+normrnd(0,5.29,[op_rent_num_i,rb_up_num,tp_num])-20)/10);
        %         case 250
        %             op_gain_up = (1/10).^((83.21+20.09*log10(Tp_bs_distance)+normrnd(0,8.41,[op_rent_num_i,rb_up_num,tp_num])-20)/10);
        %         case 270
        %             op_gain_up = (1/10).^((83.24+18.21*log10(Tp_bs_distance)+normrnd(0,9.61,[op_rent_num_i,rb_up_num,tp_num])-20)/10);
        %         case 290
        %             op_gain_up = (1/10).^((85.35+18.31*log10(Tp_bs_distance)+normrnd(0,10.24,[op_rent_num_i,rb_up_num,tp_num])-20)/10);
        %         case 310
        %             op_gain_up = (1/10).^((80.56+17.51*log10(Tp_bs_distance)+normrnd(0,7.84,[op_rent_num_i,rb_up_num,tp_num])-20)/10);
        %         case 330
        %             op_gain_up = (1/10).^((85.02+17.45*log10(Tp_bs_distance)+normrnd(0,9,[op_rent_num_i,rb_up_num,tp_num])-20)/10);
        %     end
        % tp_bs_v_gain = (reshape(total_power.*mean(op_gain_down,2),op_rent_num_i,tp_num))';
        
        %total_power_ = repmat(total_power,1,1,tp_num);
        %aa = total_power_(Tp_bs_distance_index);
        %aaaa = total_power(Tp_bs_distance_index).*mean(op_gain_down,2);
        %tp_bs_v_gain = (reshape(total_power.*mean(op_gain_down,2),s_num,tp_num))';
        %aaa=total_power_(Tp_bs_distance_index).*mean(op_gain_down,2);
%         tmp = zeros(s_num,1,tp_num);
%         for j=1:tp_num
%             tmp(:,:,j) = total_power_(Tp_bs_distance_index(:,:,j),:,j);
%         end
%         bbb = total_power_(Tp_bs_distance_index);
       
        [~,alo] = max(total_power(Tp_bs_distance_index).*mean(op_gain_down,2),[],1);
        %aaa = Tp_bs_distance_index(alo);
        %aaaa=Tp_bs_distance_index(aaa);
        %gain_ = (reshape(gain,1,tp_num))';
        %cc = Tp_bs_distance_index(alo);
       
%         tp_bs_allo_i = zeros(tp_num,1);
%         for j = 1:tp_num
%             tp_bs_allo_i(j) = Tp_bs_distance_index(alo(1,1,j),:,j);
%         end
 
        Idx = cat(1,Tp_bs_distance_index(:));
        tp_bs_allo_i=Idx(cat(1,alo(:))+[0:s_num:(tp_num-1)*s_num]');
 
        %tp_bs_allo_i = (reshape(Tp_bs_distance_index(alo),1,tp_num))';
        %tp_bs_v_gain = (reshape(total_power_(Tp_bs_distance_index).*mean(op_gain_down,2),s_num,tp_num))';
        %[~,tp_bs_allo_i] = max(tp_bs_v_gain,[],2);    
        
        %      end
        
        %%
        tps_bs_allo(row_tp) = tp_bs_allo_i;
        tps_rate_i_down = zeros(tp_num,1);
        tps_rate_i_up = zeros(tp_num,1);
         
        for j = 1:1:op_rent_num_i
            % [total_power,~,~] = ac_decode(op_rent_bs_ac_i(j,:));
            p = total_power(j)/rb_down_num;
            tps_j = find(tp_bs_allo_i == j);
            
            if isempty(tps_j)
                continue
            end
            [n,~]=size(tps_j);
            tps_j_gain_down = zeros(n,rb_down_num);
            tps_j_gain_up = zeros(n,rb_up_num);
            for k = 1:1:n
                tmp = find(Tp_bs_distance_index(:,:,tps_j(k,1)) == j);
                tps_j_gain_down(k,:) = op_gain_down(tmp,:,tps_j(k,1));
                tps_j_gain_up(k,:) = op_gain_up(tmp,:,tps_j(k,1));
            end
            %tps_j_num = length(tps_j);
            %aa = find(Tp_bs_distance_index(:,:,tps_j)==j);
            %tps_j_gain_down =  op_gain_down(aaaa,:,tps_j);
            %bbb=rem(rem((aa-1),5*1),5)+1;
            %tps_j_gain_up = op_gain_up(aaaa,:,tps_j);
            %tps_j_num = length(tps_j);
            %tps_j_gain_down = (reshape(op_gain_down(j,:,tps_j),rb_down_num,tps_j_num))';
            %tps_j_gain_up   = (reshape(op_gain_up(j,:,tps_j),rb_up_num,tps_j_num))';
  

            for m = 1:1:rb_down_num
                allo_tp = [0,0];
                allo_tp_g = [0,0];
                tmp_gain = tps_j_gain_down(:,m);
                if all(tps_rate_i_down(tps_j)>=rate_down_min) 
                   
                    tmp_gain = tps_j_gain_down(:,m);
                    for t = 1:1:2
                        [~,index] = max(tmp_gain); 
                        allo_tp_g(t) = tmp_gain(index);
                        tmp_gain(index) = -tmp_gain(index); 
                        tp_index = tps_j(index);
                        allo_tp(t) = tp_index;
                        if length(tmp_gain) == 1
                            break
                        end
                    end
                else
                   
                    exceed = find(tps_rate_i_down(tps_j)>=rate_down_min); 
                    tmp_gain(exceed) = -tmp_gain(exceed); 
                    for t = 1:1:2
                        if all(tmp_gain < 0) 
                            break
                        end
                        [~,index] = max(tmp_gain); 
                        allo_tp_g(t) = tmp_gain(index);
                        tmp_gain(index) = -tmp_gain(index); 
                        tp_index = tps_j(index);
                        allo_tp(t) = tp_index;
                    end
                end
                
                if length(find(allo_tp>0)) == 2 

                    tp1_power = p*(allo_tp_g(2)/sum(allo_tp_g)); 
                    tp2_power = p*(allo_tp_g(1)/sum(allo_tp_g)); 
                    SINR_1 = ((allo_tp_g(1))^2*tp1_power)/(1e-26); 
                    SINR_2 = ((allo_tp_g(2))^2*tp2_power)/((allo_tp_g(2))^2*tp1_power+1e-26);
                    r1 = 180000*log2(1+SINR_1);
                    r2 = 180000*log2(1+SINR_2);
                    tps_rate_i_down(allo_tp(1)) = tps_rate_i_down(allo_tp(1))+r1; 
                    tps_rate_i_down(allo_tp(2)) = tps_rate_i_down(allo_tp(2))+r2; 
                else
                   
                    tp1_power = p; 
                    SINR_1 = ((allo_tp_g(1))^2*tp1_power)/(1e-26); 
                    r1 = 180000*log2(1+SINR_1);
                    tps_rate_i_down(allo_tp(1)) = tps_rate_i_down(allo_tp(1))+r1;
                end
            end
            
            for m = 1:1:rb_up_num
                allo_tp = [0,0];
                allo_tp_g = [0,0];
                tmp_gain = tps_j_gain_up(:,m);
                p_up = 1;
                if all(tps_rate_i_up(tps_j)>=rate_up_min)
                    tmp_gain = tps_j_gain_up(:,m);
                    for t = 1:1:2
                        [~,index] = max(tmp_gain);
                        allo_tp_g(t) = tmp_gain(index);
                        tmp_gain(index) = -tmp_gain(index);
                        tp_index = tps_j(index);
                        allo_tp(t) = tp_index;
                        if length(tmp_gain) == 1
                            break
                        end
                    end
                else
                    exceed = find(tps_rate_i_up(tps_j)>=rate_up_min);
                    tmp_gain(exceed) = -tmp_gain(exceed);
                    for t = 1:1:2
                        if all(tmp_gain < 0)
                            break
                        end
                        [~,index] = max(tmp_gain);
                        allo_tp_g(t) = tmp_gain(index);
                        tmp_gain(index) = -tmp_gain(index);
                        tp_index = tps_j(index);
                        allo_tp(t) = tp_index;
                    end
                    
                end
                
                if length(find(allo_tp>0)) == 2
                    SINR_1 = ((allo_tp_g(1))^2*p_up)/(1e-26);
                    SINR_2 = ((allo_tp_g(2))^2*p_up)/((allo_tp_g(1))^2*p_up+1e-26);
                    r1 = 180000*log2(1+SINR_1);
                    r2 = 180000*log2(1+SINR_2);
                    tps_rate_i_up(allo_tp(1)) = tps_rate_i_up(allo_tp(1))+r1;
                    tps_rate_i_up(allo_tp(2)) = tps_rate_i_up(allo_tp(2))+r2;
                else
                    SINR_1 = ((allo_tp_g(1))^2*p_up)/(1e-26);
                    r1 = 180000*log2(1+SINR_1);
                    tps_rate_i_up(allo_tp(1)) = tps_rate_i_up(allo_tp(1))+r1;
                end
            end
        end
        tps_down_rate(row_tp) = tps_rate_i_down; 
        tps_up_rate(row_tp) = tps_rate_i_up;
        
    end
end


end
