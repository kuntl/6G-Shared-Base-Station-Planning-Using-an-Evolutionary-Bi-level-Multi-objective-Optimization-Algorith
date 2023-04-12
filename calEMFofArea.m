function EMF = calEMFofArea(op_rent_bs,op_rent_bs_ac,parameter)

op_num = parameter.op_num;
op_rent_num = parameter.op_rent_num;
op_frequency = parameter.op_frequency;
area = parameter.area;
interval = 50;
grid_point = cell((area/interval+1),(area/interval+1));
% Grid_point = zeros((area/interval+1)*(area/interval+1),2);

E_tot = zeros((area/interval+1),(area/interval+1));
for i = 1:1:(area/interval+1)
    for j = 1:1:(area/interval+1)
        grid_point{i,j} = [(i-1)*interval,(j-1)*interval]; 
%         Grid_point((i-1)*(area/interval+1)+j,:) = [(j-1)*interval,(i-1)*interval];
    end
end
E = cell(op_num,1);
for i = 1:1:op_num
    if i >=2
        tmp_num = op_rent_num(1:(i-1),:);
        row_bs = (sum(tmp_num(:))+1):(sum(tmp_num(:))+op_rent_num(i));
    else
        row_bs = 1:op_rent_num(i);
    end
    %row_bs = ((i-1)*op_rent_num+1):((i-1)*op_rent_num+op_rent_num);
    op_rent_bs_ac_i = op_rent_bs_ac(row_bs,:);
    op_rent_bs_i = op_rent_bs(row_bs,:);
    op_frequency_i = op_frequency(i);   
    [total_power,height,~] = ac_decode(op_rent_bs_ac_i);
    total_power = reshape(total_power,1,1,op_rent_num(i));
    height = reshape(height,1,1,op_rent_num(i));
    Op_rent_bs_i = reshape( op_rent_bs_i',1,2,op_rent_num(i));
    EIRP = 10*log(total_power)+10; 
    
    tmp = Op_rent_bs_i -cat(1,grid_point{:});
    op_bs_point_distance = sqrt((height-1.5).^2+sum(tmp.^2,2));
    op_bs_point_distance = reshape(op_bs_point_distance,[size(grid_point),op_rent_num(i)]);
    
    if gpuDeviceCount>0
          op_bs_point_distance = gpuArray(single(op_bs_point_distance));
          switch op_frequency_i
            case 230
            op_bs_point_pl = 81.10+21.61*log10(op_bs_point_distance)+5.29*gpuArray.randn([(area/interval+1),(area/interval+1),op_rent_num(i)])-10;
            case 250
            op_bs_point_pl = 83.21+20.09*log10(op_bs_point_distance)+8.41*gpuArray.randn([(area/interval+1),(area/interval+1),op_rent_num(i)])-10;
            case 270
            op_bs_point_pl = 83.24+18.21*log10(op_bs_point_distance)+9.61*gpuArray.randn([(area/interval+1),(area/interval+1),op_rent_num(i)])-10;
            case 290
            op_bs_point_pl = 85.35+18.31*log10(op_bs_point_distance)+10.24*gpuArray.randn([(area/interval+1),(area/interval+1),op_rent_num(i)])-10;
            case 310
            op_bs_point_pl = 80.56+17.51*log10(op_bs_point_distance)+7.84*gpuArray.randn([(area/interval+1),(area/interval+1),op_rent_num(i)])-10;
            case 330
            op_bs_point_pl = 85.02+17.45*log10(op_bs_point_distance)+9*gpuArray.randn([(area/interval+1),(area/interval+1),op_rent_num(i)])-10;
          end
        op_bs_point_pl = gather(op_bs_point_pl);
    else
    switch op_frequency_i
            case 230
            op_bs_point_pl = 81.10+21.61*log10(op_bs_point_distance)+normrnd(0,5.29,[(area/interval+1),(area/interval+1),op_rent_num(i)])-10;
            case 250
            op_bs_point_pl = 83.21+20.09*log10(op_bs_point_distance)+normrnd(0,8.41,[(area/interval+1),(area/interval+1),op_rent_num(i)])-10;
            case 270
            op_bs_point_pl = 83.24+18.21*log10(op_bs_point_distance)+normrnd(0,9.61,[(area/interval+1),(area/interval+1),op_rent_num(i)])-10;
            case 290
            op_bs_point_pl = 85.35+18.31*log10(op_bs_point_distance)+normrnd(0,10.24,[(area/interval+1),(area/interval+1),op_rent_num(i)])-10;
            case 310
            op_bs_point_pl = 80.56+17.51*log10(op_bs_point_distance)+normrnd(0,7.84,[(area/interval+1),(area/interval+1),op_rent_num(i)])-10;
            case 330
            op_bs_point_pl = 85.02+17.45*log10(op_bs_point_distance)+normrnd(0,9,[(area/interval+1),(area/interval+1),op_rent_num(i)])-10;
    end
     end
    E{i}= 10.^((EIRP-43.15+20*log(op_frequency_i)-op_bs_point_pl)/20); 
end

for i = 1:1:op_num  
    E_tot = E_tot +sum(E{i}.^2,3);
end
E_tot = sqrt(E_tot);
E_tot = sort(E_tot(:));
E50_index = round(length(E_tot)*0.5);
E95_index = round(length(E_tot)*0.95);
EMF = 0.5*E_tot(E50_index)+0.5*E_tot(E95_index);
end