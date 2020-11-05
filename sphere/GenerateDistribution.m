
function D = GenerateDistribution(N)
sigma = 0.05;
myu = 0;
x = linspace(0, 1, N );
D = NormalDistribution(x, myu, sigma);
D = D / max( D );
end

function D = NormalDistribution(x, myu, sigma)
D = 1./sigma./sqrt(2*pi)*exp(-(x-myu).^2/2./sigma.^2);
end
