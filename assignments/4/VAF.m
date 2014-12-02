function [vaf]=VAF(u,usim)

            %Time Domain Validation of fits: VAFs 
            vaf=( 1  - sum ((abs(u(1:end)-usim(1:end))).^2)./ sum( (abs(u(1:end))).^2)  )*100;
