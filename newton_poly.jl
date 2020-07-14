using PyPlot
k = 51
#Chebyshev nodes
xp = Vector(undef,k)
x = Vector(undef,k)
for i = 1:k
    xp[i] = cos(π*(2i-1)/(2(k+1)))
    x[i] = cos(π*(i/k))
end
# Equispaced nodes
#xp = Vector(LinRange(-1,1,k))
f(x) = 1/sqrt(1-abs(x))
w(s) = 1/sqrt(1-s^2)
#A = zeros(k,k)
# for i=2:k
#      for j = 1:k #-i+1
#          #dd[j,i] = (dd[j+1,i-1]-dd[j,i-1])/(xp[j+i-1]-xp[j])
#          A[i,j] = 1/sqrt(1-(xp[i]-x[j])^2)
#      end
# end
W = π*w.(xp)./k
y = zeros(k)
#A = dd[1,:]'
#for i = 1:k
    for j = 1:k
        y[j] += W[j]*f(xp[j])#/(x[i]-xp[j])
    end
#end
println(cumsum(y,dims=1))
#x = LinRange(-1,1,m)
# x1 = Vector(ones(k))
# y = A[1].*x1
# for i = 2:k
#     for j = 1:k
#         x1[j] /=  (x[j]-xp[i-1])
#         y[j] += A[i]*x1[j]
#     end
# end
#plot(x,f.(x))
#plot(x,w.(xp))
#ylim(0,1)
