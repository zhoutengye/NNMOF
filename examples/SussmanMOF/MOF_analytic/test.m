V = 1/3;

dd= 0.01;
x = linspace(dd,1,100);

X = [];
Y = [];
Z = [];


for i = 1:100
	x1 = x(i);
	V2 = V/x1;
	y = linspace(V2,1,100);
	for j = 1:100
		y1 = y(j);
		z1 = V/x1/y1;
		X = [X,x1/3];
		Y = [Y,y1/3];
		Z = [Z,z1/3];
	end
	
end

plot3(X,Y,Z)
axis([0,1/3,0,1/3,0,1/3])