% smooth TV for Bernoulli case

function[d] = smooth_tv_b(x,y,px,py,n,q,L)

% x, y are the location with positive point mass
% px, py are the point mass at the locations in x, y
% n is the dimension of the matrix
% q is the probability parameter in Bernoulli (true parameter or estimated q)
% L is the smooth parameter

t_seq = -1/2:1/L:1/2;
u_seq = t_seq.*sqrt(n*q*(1-q)) + n*q;

fx = zeros(L,1);
fy = zeros(L,1);

for l = 2:1:(L+1)

    if u_seq(l) >= min(x)
        if u_seq(l-1) >= min(x)
            fx(l-1) = sum(px( logical( x <= u_seq(l) ) ) ) - sum(px( logical( x <= u_seq(l-1) ) ) );
        else 
            fx(l-1) = sum(px( logical( x <= u_seq(l) ) ) );
        end
    end

    if u_seq(l) >= min(y)
        if u_seq(l-1) >= min(y)
            fy(l-1) = sum(py( logical( y <= u_seq(l) ) ) ) - sum(py( logical( y <= u_seq(l-1) ) ) );
        else 
            fy(l-1) = sum(py( logical( y <= u_seq(l) ) ) );
        end
    end
end

d = sum(abs(fx - fy));

