function [total_power,height,cost]=ac_decode(config)
% total_power_type = 500*[1,2,3,4,5,6,7,8,9,10];
% total_power_type_cost = [1,2,3,4,5,6,7,8,9,10];
% height_type = [5,10,15,20];
% height_type_cost = [1,2,3,4];
% height = (height_type(config(:,1)))';
% total_power = (total_power_type(config(:,2)))';
% cost = (total_power_type_cost(config(:,2)) + height_type_cost(config(:,1)))';

height = config(:,1,:)*5;
total_power = config(:,2,:)*500;
cost = config(:,2,:) + config(:,1,:);

end