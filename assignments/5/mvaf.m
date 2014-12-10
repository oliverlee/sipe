function [vaft]=mvaf(y,y0,N)
dummy=ceil(N/2);
vaft=zeros(length(y),1);
for jj=dummy+1:length(y)-dummy;
    vaft(jj)=vaf(y(jj-dummy:jj+dummy),y0(jj-dummy:jj+dummy));
end

