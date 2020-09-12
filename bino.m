% So = 40;
% K = 45;
% r = 0.08;
% T = 5/12;
% sigma = 0.30;
% N = 25:250;
price = zeros(size(N));
K = 100; r = 0.1; sigma = 0.30;
T = 1;
N =101;
M = 100;
y = linspace(0,6,N);
So = exp(y);
for k = 1:N
 delta = T/M;
 [StockPrice, OptionPrice] = binprice(So, K, r, T, delta, sigma, 0);
 price(k) = OptionPrice(1,1);
end
plot(N, price);
xlabel('Iterations');
ylabel('Option Value');
title('Binomial Pricing');