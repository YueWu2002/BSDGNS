function [val] = minmod_3(a,b1,b2)
% minmod function with 3 inputs
% vectorization supported

a = a(:);
b1 = b1(:);
b2 = b2(:);

val = zeros(size(a));
flag1 = a > 0.0;
flag2 = a < 0.0;
val(flag1) = max(min([a(flag1), b1(flag1), b2(flag1)], [], 2), 0.0);
val(flag2) = min(max([a(flag2), b1(flag2), b2(flag2)], [], 2), 0.0);
end