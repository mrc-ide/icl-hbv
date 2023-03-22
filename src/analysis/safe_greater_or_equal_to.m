function out = safe_greater_or_equal_to(a,b)
% Tests to see if a>=b without errors causing by floating points exceptions

out = or(logical(a>b), logical(abs(a-b)<0.00001));