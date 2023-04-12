classdef SharedBS6G_2000 < PROBLEM
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
        function obj = SharedBS6G_2000()
            fullpath = mfilename('fullpath'); 
            [path,~]=fileparts(fullpath);
            try
                load([path,'\Setting_2000.mat']);
            catch

                rb_down_num = 300;
                rb_up_num = 150;
                cs_num = 500;
                area = 5000;  
                CS = area*rand(cs_num,2); 
                r_tmp = rand(cs_num,1);
                CS_COST = 30+50*r_tmp;
                Bh_C_u = 0.00000001*r_tmp;  
                op_frequency = [230,250,270,290,310,330];
                CS_RENT = 1.5*CS_COST; 
                op_num = 3;
                tp_num = 2000;
                TP =area*rand(op_num*tp_num,2);
                rate_up_min = 1e4;
                rate_down_min = 1e5;     
                height_type = 4;       
                total_power_type = 10;   
                phi = 0.6;
                BGmax = 40000; 
                e  = 5000;     
                Fcov_min = 0.95;
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
            bs_rent_count = sum(y_site,3);
            
            D = 1+0.2*bs_rent_count;
            CS_RENT = CS_RENT'./D; 
            Y =zeros(size(y));
            Y(y_site) = 1; 
            ops_rent_sum_i = reshape(sum(Y.*CS_RENT,2),N,op_num);  
            op_rent_num = (reshape(sum(Y,2),N,op_num))'; 
            
            Index_ops_CS = CS_Index.*Y;
            Index_ops_CS = (reshape(Index_ops_CS,N,cs_num*op_num))';                          
            
            Configure = zeros(cs_num*op_num,2,N); 
            Heights = ceil(lDecs/total_power_type);
            Configure(:,1,:) = Heights';    
            Configure(:,2,:) = (lDecs - max(0,Heights-1)*total_power_type)'; 
            
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
                       Parameter.Bh_BS_u = Bh_BS_u;           
                       
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
                    fcov_i = reshape((sum(cov_i)/tp_num),N,op_num);
                    fcap_i = reshape(sum(TPS_up_rateS)+sum(TPS_up_rateS),N,op_num);
                    cost = reshape(cost,cs_num,op_num,N);
                    AC_sum_i = sum(cost);
                    AC_sum_i = reshape(AC_sum_i(:),op_num,N)';
                    
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
                       
                       [~,Parameter.index_ops_bs] = ismember(index_ops_cs,unique(index_ops_cs)); 

                       op_rent_bs_ac = Configure(:,:,i);
                       op_rent_bs_ac = op_rent_bs_ac(site,:);
                       
                       EMF(i) = calEMFofArea(op_rent_bs,op_rent_bs_ac,Parameter);
                       
                       
                       Bh_BS_u = Bh_C_u(x(i,:)>0);
                       Parameter.Bh_BS_u = Bh_BS_u;                    
                     
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
                   fcov_i = reshape((sum(cov_i)/tp_num),N,op_num); 
                   fcap_i = reshape(sum(TPS_up_rateS)+sum(TPS_up_rateS),N,op_num);
                   cost = reshape(cost,cs_num,op_num,N);
                   AC_sum_i = sum(cost); 
                   AC_sum_i = reshape(AC_sum_i(:),op_num,N)';
                    
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
                    phi = obj.Parameter.phi;
                    BGmax = obj.Parameter.BGmax; 
                    e  = obj.Parameter.e; 
                    cs_num = obj.Parameter.cs_num;
                    op_num = obj.Parameter.op_num;
                    EMF = cat(1,parameters.EMF);
                    InC = cat(1,parameters.InC);
                    BhC = cat(1,parameters.BhC);
                    u_Con(:,1) = phi - sum(uDecs,2)/size(uDecs,2); 
                    u_Con(:,2) = InC + BhC - BGmax; 
                    y = reshape(y,N,cs_num,op_num);
                    y = y.*uDecs;                    
                    u_Con(:,3) = sum(x,2) - sum(sum(y,3)>0,2);
                    u_Con(:,4) = EMF-e; 
                    l_Con = [];
                case 'lower'
                    u_Con = [];
                    Fcov_min = obj.Parameter.Fcov_min;
                    fcov = cat(1,parameters.fcov);
                    l_Con = Fcov_min - fcov; 
                case 'bilevel'
                    % upper-level
                    phi = obj.Parameter.phi;
                    BGmax = obj.Parameter.BGmax;
                    e  = obj.Parameter.e;    
                    cs_num = obj.Parameter.cs_num;
                    op_num = obj.Parameter.op_num;
                    EMF = cat(1,parameters.EMF);
                    InC = cat(1,parameters.InC);
                    BhC = cat(1,parameters.BhC);
                    u_Con(:,1) = phi - sum(uDecs,2)/size(uDecs,2); 
                    u_Con(:,2) = InC + BhC - BGmax; 
                    y = reshape(y,N,cs_num,op_num);
                    y = y.*uDecs;                    
                    u_Con(:,3) = sum(sum(y,3)>0,2)-sum(x,2);
                    u_Con(:,4) = EMF-e;
                    
                    %lower-level
                    Fcov_min = obj.Parameter.Fcov_min;
                    fcov = cat(1,parameters.fcov);
                    l_Con = Fcov_min - fcov; 
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