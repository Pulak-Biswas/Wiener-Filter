function [Dac,lags,otpshot]=WSDF(inpshot,Flen,WN,ts)
% %-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% % WIENER SPIKING DECONVOLUTION FILTER -%-  PULAK BISWAS   %
% % Example: [Dac,lags,otpshot]=WSDF(inpshot,Flen,WN,ts)    %
% %       INPUTS     inpshot: seismic shot gather(s)        %
% %                  Flen: Filter length                    %
% %                  WN: White noise percentage             %
% %                  ts: time sampling interval             %
% %      OUTPUTS     Dac: Autocorrelation fnc               %
% %                  lags: Autocorrelation lags             %
% %                  otpshot: actual output                 %
% %-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

box=size(inpshot);               % size of input data
Dwav=[1,zeros(1,box(1)-1)];      % desired wavelet

Dac=zeros(Flen,box(2));          
Dcc=zeros(1,Flen);
for i=1:box(2)
    corr=zeros(1,Flen);
    for j=1:Flen
        for k=1:box(1)
            if k+j <= box(1)+1
                corr(j)=corr(j)+inpshot(k,i)*inpshot(k+j-1,i);
            end
        end
    end
    Dac(:,i)=corr; 
    corr=zeros(1,Flen);
    for j=1:Flen
        for k=1:box(1)
            if k+j <= box(1)+1
                corr(j)=corr(j)+inpshot(k,i)*Dwav(k+j-1)';
            end
        end
    end
    Dcc=Dcc+corr;
end
lags=0:ts:(Flen-1)*ts;

Dcc=abs(Dcc);
Dac2=sum(Dac,2);
Dac2(1)=Dac2(1)+WN/100;
Tpltz=toeplitz(Dac2);

otpshot=zeros(box(1)+Flen-1,box(2));
hfilt=Dcc/(Tpltz);
for i=1:box(2)
    otpshot(:,i)=conv(inpshot(:,i),hfilt);
end

otpshot=otpshot(1:box(1),:);


