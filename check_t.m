function tt = check_t(t)
% handle GNSS time overflow
half_week = 302400;
tt = t;
if t >  half_week, tt = t-604800; end
if t < -half_week, tt = t+604800; end
