function [week,sow] = toGPStime(jday)
%  convert from Julian Day (jday) number to GPS week and Seconds of Week

a = floor(jday+0.5);
b = a+1537;
c = floor((b-122.1)/365.25);
e = floor(365.25*c);
f = floor((b-e)/30.6001);
d = b-e-floor(30.6001*f)+rem(jday+0.5,1);
dow = rem(floor(jday + 0.5),7);%day of week
week = floor((jday-2444244.5)/7);
%second of week
sow = (rem(d,1)+dow+1)*86400;
