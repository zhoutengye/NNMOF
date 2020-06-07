clear
clc

n = 101;

x = linspace(0,1,n);
y = linspace(0,1,n);
V=zeros(n,n);
flag=zeros(n,n);

for i = 1:1:51
	for j = 1:1:51
		if(x(i)<=1/3 && y(j)<=1/3)
			V(i,j) = 9/2*x(i)*y(j);
			flag(i,j) = 1;
		else
			if(x(i)>=1/3 && x(i)<=1/2 && y(j)<=1/2)
				V(i,j) = y(j) / (0.5+6*(0.5-x(i))^2);
				flag(i,j) = 1;
				if(y(j)>V(i,j)/2+1/24/V(i,j)*(2*V(i,j)-2)^2 )
					V(i,j) = 0;
					flag(i,j) = 0;
				end
			end
			if(y(j)>=1/3 && y(j)<=1/2 && x(i)<=1/2 && flag(i,j) ==0)
				V(i,j) = x(i) / (0.5+6*(0.5-y(j))^2);
				flag(i,j) = 1;
				if(x(i)>V(i,j)/2+1/24/V(i,j)*(2*V(i,j)-2)^2 )
					V(i,j) = 0;
					flag(i,j) = 0;
				end
			end
		end

	end
end

for i = 1:1:51
	for j = 1:1:51
		if(flag(i,j)==0)
			a1 = 1;
			a2 = -3 + 4.5*(1-x(i))*(1-y(j));
			a3 = 3 - 2.25*(1-x(i)) - 2.25*(1-y(j));
			a4 = 1/8;
			r = roots([a1,a2,a3,a4]);
			for kk = 1:length(r)
				if(r(kk)>=1/2 && r(kk)<=1) 
					V(i,j) = r(kk);
				end
			end
		end
	end
end

for i = 1:1:51
	for j = 1:1:51
		V(102-i,102-j) = V(i,j);
		V(102-i,j)     = V(i,j);
		V(i,102-j)     = V(i,j);
	end
end

[X,Y] = meshgrid(x,y);

surfc(X,Y,V,'EdgeColor','none','FaceAlpha',0.5)
h_contours = findall(gcf,'Type','patch');
set(h_contours, 'LineWidth', 5)