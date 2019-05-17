sol = zeros(7,2);
w = [0,0.5,1,2,4,8,16];


for i=1:7
t = 0:80;
q_init = [0;0];
[t,q] = ode45(@ode,t,q_init);

function dqdt = ode(t,q)
q1 = 0;

q2 = 0;

dqdt_1 = 0;

dqdt_2 = 0;

dqdt = [0;0];

dqdt_1 = q2;

dqdt_2 = - 4*q1 - 5*q2 + 10*cos(w(i)*t);

dqdt = [dqdt_1; dqdt_2];

end

sol(i) = [t,y];
end


for i = 1:7

plot(sol(i));

title('w = '%d' ',w(i))

xlable('t')

ylable('q')

end