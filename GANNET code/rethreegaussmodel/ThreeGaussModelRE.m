function F = ThreeGaussModel(x,freq)
A=0.058;
F = x(1)*(...
       x(6)*exp(x(2)*(freq-x(3)).*(freq-x(3)))+exp(x(2)*(freq-x(3)+A).*(freq-x(3)+A))+exp(x(2)*(freq-x(3)-A).*(freq-x(3)-A))...
    )...
    +x(4)*(freq-x(3))+x(5);

end