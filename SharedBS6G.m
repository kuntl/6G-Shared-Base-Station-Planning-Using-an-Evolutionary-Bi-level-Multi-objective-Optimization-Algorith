classdef SharedBS6G < PROBLEM
% <problem> <DS>
% Benchmark multi-objective bi-level problem

%------------------------------- Reference --------------------------------
% Deb K , Sinha A . An Evolutionary Approach for Bilevel Multi-objective 
% Problems[J]. Communications in Computer & Information Science, 2009,
% 35:17-24.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        %% Initialization
        function obj = SharedBS6G()
            fullpath = mfilename('fullpath'); 
            [path,~]=fileparts(fullpath);
            try
                load([path,'\Setting.mat']);
            catch
               %上下行信道数量
                rb_down_num = 300;
                rb_up_num = 100;
                %随机生成cs_num个候选站址
                cs_num = 500;
                area = 3000;  %单位：米
                CS = area*rand(cs_num,2);  %候选站址
                %在每个候选站址建设基站的费用
                r_tmp = rand(cs_num,1);
                CS_COST = 30+50*r_tmp;

                %候选站址回程链路的单位费用
                Bh_C_u = 0.00000001*r_tmp;
                
                op_frequency = [230,250,270,290,310,330];

                %每个候选站址的租金
                CS_RENT = 1.5*CS_COST; %未考虑共享折扣
                
                %每个运营商随机生成tp_num个测试点
                op_num = 4;
                tp_num = 800;
                TP =area*rand(op_num*tp_num,2);
                
                %用户的最低上下行速率
                rate_up_min = 1e6;
                rate_down_min = 10e6;     
                
                height_type = 4;       %天线高度类型数量
                total_power_type = 10; %天线功率类型数量
                
                phi = 0.6;%最小基站密度
                BGmax = 4000; %预算
                e  = 1800; %电磁场暴露阈值
                
                Fcov_min = 0.95;%最小覆盖率
                
                save([path,'\Setting.mat'],'rb_down_num','rb_up_num','cs_num','area','CS','CS_COST','Bh_C_u','op_frequency','CS_RENT','op_num','tp_num','TP','rate_up_min','rate_down_min','height_type','total_power_type','phi','BGmax','e','Fcov_min');
            end
                        
            obj.Parameter.rb_down_num = rb_down_num;
            obj.Parameter.rb_up_num = rb_up_num;
            obj.Parameter.cs_num = cs_num;
            obj.Parameter.area = area;
            obj.Parameter.CS = CS;
            obj.Parameter.CS_COST = CS_COST;
            obj.Parameter.Bh_C_u = Bh_C_u;
            obj.Parameter.op_frequency = op_frequency;
            obj.Parameter.CS_RENT = CS_RENT;
            obj.Parameter.op_num = op_num;
            obj.Parameter.tp_num = tp_num;
            obj.Parameter.TP = TP;
            obj.Parameter.rate_up_min = rate_up_min;
            obj.Parameter.rate_down_min = rate_down_min;
            obj.Parameter.height_type = height_type;
            obj.Parameter.total_power_type= total_power_type;
            obj.Parameter.phi = phi;
            obj.Parameter.BGmax = BGmax;
            obj.Parameter.e = e;
            obj.Parameter.Fcov_min = Fcov_min;
            
            
            obj.Global.u_M = 2;
            obj.Global.l_M = 3;
            obj.Global.u_D = cs_num;
            obj.Global.l_D = cs_num*op_num;
            
            obj.Global.u_lower    = zeros(1,obj.Global.u_D);
            obj.Global.u_upper    = ones(1,obj.Global.u_D);
            obj.Global.l_lower    = zeros(1,obj.Global.l_D);
            obj.Global.l_upper    = height_type*total_power_type*ones(1,obj.Global.l_D);
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function [u_Obj,l_Obj,parameters] = CalObj(obj,uDecs,lDecs,parameters,type)
 
             %Loading parameters 
            Parameter = obj.Parameter;
            cs_num = Parameter.cs_num;
            op_num = Parameter.op_num;
            CS_RENT = Parameter.CS_RENT;
            total_power_type = Parameter.total_power_type;
            CS = Parameter.CS;
            tp_num = Parameter.tp_num;
            TP = Parameter.TP;
            CS_COST = Parameter.CS_COST;
            Bh_C_u = Parameter.Bh_C_u;
            
            %Decoding
            x = uDecs;
            y = lDecs;
            
            N = size(x,1);
                     
            CS_Index = repmat(1:cs_num,N,1);
            
            bs_num = sum(x,2);
            
            y = reshape(y,N,cs_num,op_num);
            y = y.*x;
            lDecs = reshape(y,N,cs_num*op_num);
            
            y_site = y>0;  
            bs_rent_count = sum(y_site,3); %每个基站上的运营商数量
            
            D = 1+0.2*bs_rent_count; %运营商在每个基站上享有的租金折扣
            CS_RENT = CS_RENT'./D;  %运营商在每个基站上的租金成本
            Y =zeros(size(y));
            Y(y_site) = 1; 
            ops_rent_sum_i = reshape(sum(Y.*CS_RENT,2),N,op_num);   %各个运营商的总租金成本
            op_rent_num = (reshape(sum(Y,2),N,op_num))';  %每个运营商租用的基站数量
            
            Index_ops_CS = CS_Index.*Y;
            Index_ops_CS = (reshape(Index_ops_CS,N,cs_num*op_num))';%运营商租用基站在候选基站中的编号                                          
            
            Configure = zeros(cs_num*op_num,2,N); %每个运营商在每个基站上的天线配置
            Heights = ceil(lDecs/total_power_type);
            Configure(:,1,:) = Heights';     %天线高度类型
            Configure(:,2,:) = (lDecs - max(0,Heights-1)*total_power_type)'; %天线功率类型
            
            [total_power,height,cost]=ac_decode(Configure);
                                         
                                
            
            if isempty(parameters)
                for i = 1:N
                    Parameter.bs_num = bs_num(i);
                    Parameter.op_rent_num = op_rent_num(:,i);
                    
                    index_ops_cs = Index_ops_CS(:,i);
                    site = index_ops_cs>0;
                    
                    index_ops_cs = index_ops_cs(site);
                    op_rent_bs = CS(index_ops_cs,:);
                    
                    op_rent_bs_ac = Configure(:,:,i);
                    op_rent_bs_ac = op_rent_bs_ac(site,:);
                    
%                     [parameters(i).TPS_down_rate,parameters(i).TPS_up_rate,parameters(i).TP_bs_allo] = calRateofTps(TP,op_rent_bs,op_rent_bs_ac,Parameter);
                    [parameters(i).TPS_down_rate,parameters(i).TPS_up_rate,parameters(i).TP_bs_allo] = calRateofTps_v3(TP,op_rent_bs,op_rent_bs_ac,Parameter);
                end               
            end
            TPS_down_rateS = cat(2,parameters.TPS_down_rate);
            TPS_up_rateS   = cat(2,parameters.TPS_up_rate);
            TP_bs_alloS    = cat(2,parameters.TP_bs_allo);
            
            switch type
                case 'upper'
                   InC = sum(repmat(CS_COST',N,1).*x,2);
                   BhC = zeros(N,1);
                   P_tot = sum(total_power);
                   EMF =zeros(N,1);
                       
                   for i = 1:N 
                       Parameter.bs_num = bs_num(i);
                       Parameter.op_rent_num = op_rent_num(:,i);
                       index_ops_cs = Index_ops_CS(:,i);
                       site = index_ops_cs>0;
                       index_ops_cs = index_ops_cs(site);                                             
                       op_rent_bs = CS(index_ops_cs,:);
                       
                       [~,Parameter.index_ops_bs] = ismember(index_ops_cs,unique(index_ops_cs)); %每个运营商租用基站的编号

                       op_rent_bs_ac = Configure(:,:,i);
                       op_rent_bs_ac = op_rent_bs_ac(site,:);
                       
                       EMF(i) = calEMFofArea(op_rent_bs,op_rent_bs_ac,Parameter);                       
                       
                       Bh_BS_u = Bh_C_u(x(i,:)>0);
                       Parameter.Bh_BS_u = Bh_BS_u; %每个基站单位回程链路建设成本                       
                       
                       BhC(i) = sum(calBhCofBSs(TPS_down_rateS(:,i),TPS_up_rateS(:,i),TP_bs_alloS(:,i),Parameter));  
                       
                       parameters(i).EMF = EMF(i);
                       parameters(i).BhC=BhC(i);
                       parameters(i).InC = InC(i);
                      
                   end
                   
                   
                   u_Obj(:,1) =  InC + BhC - sum(ops_rent_sum_i,2);
                   u_Obj(:,2) =  P_tot(:);
                   
                   
                   l_Obj = [];
                case 'lower'
                    u_Obj = [];
                    
                    rate_down_min = obj.Parameter.rate_down_min;
                    rate_up_min = obj.Parameter.rate_up_min;
                    
                    TPS_up_rateS = reshape(TPS_up_rateS,tp_num,N,op_num);
                    TPS_down_rateS = reshape(TPS_down_rateS,tp_num,N,op_num);
                    cov_i = TPS_down_rateS > rate_down_min & TPS_up_rateS > rate_up_min;
                    fcov_i = reshape((sum(cov_i)/tp_num),N,op_num); % 各个运营商的网络覆盖率
                    fcap_i = reshape(sum(TPS_up_rateS)+sum(TPS_up_rateS),N,op_num);
                    cost = reshape(cost,cs_num,op_num,N);
                    AC_sum_i = sum(cost); %各个运营商的天线总花费
                    AC_sum_i = reshape(AC_sum_i(:),op_num,N)';%运营商的天线建设成本
                    
                    fcov = min(fcov_i,[],2);
                    l_Obj(:,1) = -fcov;
                    l_Obj(:,2) = -min(fcap_i,[],2);
                    l_Obj(:,3) = max(ops_rent_sum_i+AC_sum_i,[],2);
                    
                    
                    for i = 1:N
                        parameters(i).fcov = fcov;
                    end
                    
                case 'bilevel'
                    % upper-level
                    InC = sum(repmat(CS_COST',N,1).*x,2);
                    BhC = zeros(N,1);
                    P_tot = sum(total_power);
                    EMF =zeros(N,1);
                       
                    for i = 1:N 
                       Parameter.bs_num = bs_num(i);
                       Parameter.op_rent_num = op_rent_num(:,i);
                       index_ops_cs = Index_ops_CS(:,i);
                       site = index_ops_cs>0;
                       index_ops_cs = index_ops_cs(site);                                             
                       op_rent_bs = CS(index_ops_cs,:);
                       
                       [~,Parameter.index_ops_bs] = ismember(index_ops_cs,unique(index_ops_cs)); %每个运营商租用基站的编号

                       op_rent_bs_ac = Configure(:,:,i);
                       op_rent_bs_ac = op_rent_bs_ac(site,:);
                       
                       EMF(i) = calEMFofArea(op_rent_bs,op_rent_bs_ac,Parameter);
                       
                       
                       Bh_BS_u = Bh_C_u(x(i,:)>0);
                       Parameter.Bh_BS_u = Bh_BS_u; %每个基站单位回程链路建设成本                       
                     
                       BhC(i) = sum(calBhCofBSs(TPS_down_rateS(:,i),TPS_up_rateS(:,i),TP_bs_alloS(:,i),Parameter));     
                       
                       
                       parameters(i).EMF = EMF(i);
                       parameters(i).BhC=BhC(i);
                       parameters(i).InC = InC(i);
                       
                    end                  
             
                   u_Obj(:,1) =  InC + BhC - sum(ops_rent_sum_i,2);
                   u_Obj(:,2) =  P_tot(:);

                   
                   %lower-level
                   rate_down_min = obj.Parameter.rate_down_min;
                   rate_up_min = obj.Parameter.rate_up_min;
                    
                   TPS_up_rateS = reshape(TPS_up_rateS,tp_num,N,op_num);
                   TPS_down_rateS = reshape(TPS_down_rateS,tp_num,N,op_num);
                   cov_i = TPS_down_rateS > rate_down_min & TPS_up_rateS > rate_up_min;
                   fcov_i = reshape((sum(cov_i)/tp_num),N,op_num); % 各个运营商的网络覆盖率
                   fcap_i = reshape(sum(TPS_up_rateS)+sum(TPS_up_rateS),N,op_num);
                   cost = reshape(cost,cs_num,op_num,N);
                   AC_sum_i = sum(cost); %各个运营商的天线总花费
                   AC_sum_i = reshape(AC_sum_i(:),op_num,N)';%运营商的天线建设成本
                    
                   fcov = min(fcov_i,[],2);
                   l_Obj(:,1) = -fcov;
                   l_Obj(:,2) = -min(fcap_i,[],2);
                   l_Obj(:,3) = max(ops_rent_sum_i+AC_sum_i,[],2);
                   
                   
                   for i = 1:N
                            parameters(i).fcov = fcov;
                   end
                   
            end
        end
        %% Calculate constraint violations
        function [u_Con,l_Con] = CalCon(obj,uDecs,lDecs,parameters,type)
            x = uDecs;
            y = lDecs;
            N = size(x,1);
            
            switch type
                case 'upper'
                    phi = obj.Parameter.phi;%最小基站密度
                    BGmax = obj.Parameter.BGmax; %预算
                    e  = obj.Parameter.e; %电磁场暴露阈值
                    
                    cs_num = obj.Parameter.cs_num;
                    op_num = obj.Parameter.op_num;
                    EMF = cat(1,parameters.EMF);
                    InC = cat(1,parameters.InC);
                    BhC = cat(1,parameters.BhC);
                    
                    %基站密度不能低于phi
                    u_Con(:,1) = phi - sum(uDecs,2)/size(uDecs,2); 
                    
                    %建设总成本不能超过预算
                    u_Con(:,2) = InC + BhC - BGmax; 
                    
                    %每个基站至少被一个运营商选择
                    y = reshape(y,N,cs_num,op_num);
                    y = y.*uDecs;                    
                    u_Con(:,3) = sum(x,2) - sum(sum(y,3)>0,2);
                    
                    %服务区域内电磁场暴露不能超过阈值
                    u_Con(:,4) = EMF-e; 
                    
                    l_Con = [];
                case 'lower'
                    u_Con = [];
                    
                    Fcov_min = obj.Parameter.Fcov_min;
                    fcov = cat(1,parameters.fcov);
                    l_Con = Fcov_min - fcov; %运营商需要满足最低的覆盖率
                case 'bilevel'
                    % upper-level
                    phi = obj.Parameter.phi;%最小基站密度
                    BGmax = obj.Parameter.BGmax; %预算
                    e  = obj.Parameter.e; %电磁场暴露阈值
                    
                    cs_num = obj.Parameter.cs_num;
                    op_num = obj.Parameter.op_num;
                    EMF = cat(1,parameters.EMF);
                    InC = cat(1,parameters.InC);
                    BhC = cat(1,parameters.BhC);
                    
                    %基站密度不能低于phi
                    u_Con(:,1) = phi - sum(uDecs,2)/size(uDecs,2); 
                    
                    %建设总成本不能超过预算
                    u_Con(:,2) = InC + BhC - BGmax; 
                    
                    %每个基站至少被一个运营商选择
                    y = reshape(y,N,cs_num,op_num);
                    y = y.*uDecs;                    
                    u_Con(:,3) = sum(sum(y,3)>0,2)-sum(x,2);
                    
                    %服务区域内电磁场暴露不能超过阈值
                    u_Con(:,4) = EMF-e;
                    
                    %lower-level
                    Fcov_min = obj.Parameter.Fcov_min;
                    fcov = cat(1,parameters.fcov);
                    l_Con = Fcov_min - fcov; %运营商需要满足最低的覆盖率
            end
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
%             try
%                 load('Problems\DS\DS1.mat');
%             catch
                P = [];
               
%             end
        end
        
        function P = lower_PF(obj,y)
%             try
%                 load('Problems\DS\DS1.mat');
%             catch
                P = [];
               
%             end
        end
    end
end