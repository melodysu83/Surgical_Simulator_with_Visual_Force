function [ output_args ] = deform(tissue_param)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    persistent Ko K1 K2 ko k1 k2 b1 b2 A B;
    if isempty(Ko)
        Ko = tissue_param(1);
        K1 = tissue_param(2);
        K2 = tissue_param(3);
        ko = tissue_param(4);
        k1 = tissue_param(5);
        k2 = tissue_param(6);
        b1 = tissue_param(7);
        b2 = tissue_param(8);
        A = tissue_param(9);
        B = tissue_param(10);
    end

end

