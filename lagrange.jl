using PyPlot
f(x) = x^2
k = 4
m = 11
xp = LinRange(-1,1,k)
yp = f.(xp)
x = LinRange(-1,1,m)
w = ones(k,m)
for i = 1:k
    for j = 1:k
        if i != j
            w[i,:] = w[i,:].*(x.-xp[j])./(xp[i]-xp[j])
        end
    end
end
y = zeros(m)
for i = 1:k
    for j = 1:m
        y[j] += w[i,j]*yp[i]
    end
end
plot(x,f.(x))
plot(x,y)
