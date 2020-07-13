function [data, y] = int_u(prob, data, u) %#ok<INUSL>

ybp = u(1:end-1);
T   = u(end);
ycn = data.Wda*ybp;
y   = 0.5*T*data.wts1*ycn/data.ddaecoll.NTST;

end