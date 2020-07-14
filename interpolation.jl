using PyPlot
f(x) = 1/(1+x^2)
k = 8
m = 1000
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
