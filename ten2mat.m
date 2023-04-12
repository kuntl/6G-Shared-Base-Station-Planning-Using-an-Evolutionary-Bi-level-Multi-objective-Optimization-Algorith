function mat = ten2mat( ten,k)
% Tensor unfolding for ten with size of n1*n2*...*nd
% dim = (n1,n2,...,nd)
% k=1,2,...,d
%% Example
% ten(:,:,1) =                          ten(:,:,2) =
% 
%      1     5     9                            13    17    21
%      2     6    10                            14    18    22
%      3     7    11                            15    19    23
%      4     8    12                            16    20    24

% >> ten2mat(ten,1)
 
% ans =
 
%      1     5     9    13    17    21
%      2     6    10    14    18    22
%      3     7    11    15    19    23
%      4     8    12    16    20    24

% >> ten2mat(ten,2)

% ans =
 
%      1     2     3     4    13    14    15    16
%      5     6     7     8    17    18    19    20
%      9    10    11    12    21    22    23    24

% >> ten2mat(ten,3)

% ans =

%      1     2     3     4     5     6     7     8     9    10    11    12
%     13    14    15    16    17    18    19    20    21    22    23    24

%%
dim = size(ten);
mat = reshape(permute(ten,[k,1:k-1,k+1:length(dim)]),dim(k),[]);
end
