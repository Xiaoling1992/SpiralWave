xs= [54, 59, 66, 71, 71, 76, 79, 73, 69, 69, 56];
ys= [68, 55, 56, 60, 66, 71, 80, 75, 68, 60, 55];
plot(xs, ys)
dis= 0;
for i= 2: 1: length(xs)
    dis= dis+ sqrt( ( xs(i)-xs(i-1) )^2+ ( ys(i)-ys(i-1) )^2 )
end

speed= dis/ (4* length(xs) )