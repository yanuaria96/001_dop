function s = Integr(x,y)
s = sum(0.5 * (y(2:end) + y(1:end-1)) .* diff(x));
end