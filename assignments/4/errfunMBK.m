function [e,Hef] = errfunMBK(P, frf, fvec,mCoh)
   
    M    = P(1);
    B    = P(2);
    K    = P(3);
           
    s=2*pi*1j*fvec;
    s2=(-4*pi^2*fvec.^2);

    Hef=1./(M*s2+B.*s+K);

    e=abs(log(frf./Hef)); % 
end



