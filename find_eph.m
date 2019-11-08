function icol = find_eph(Eph,sv,time)
% Finds the proper column in ephemeris array

icol = 0;
isat = find(Eph(1,:) == sv);
n = size(isat,2);
if n == 0
   return
end
icol = isat(1);
dtmin = Eph(21,icol)-time;
% output the eph closest to the input 'time'
for t = isat
   dt = Eph(21,t)-time;
   if dt < 0
      if abs(dt) < abs(dtmin)
         icol = t;
         dtmin = dt;
      end
   end
end
