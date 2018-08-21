function [ normX ] = claculate_norm(X,K,PARFOR_FLAG)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    normX=0;
    if (PARFOR_FLAG)
        parfor k = 1:K
      
         normX = normX + sum(sum( (X{k}).^2));
        end
    else
        for k = 1:K
         normX = normX + sum(sum( (X{k}).^2));
        end
    end


end

