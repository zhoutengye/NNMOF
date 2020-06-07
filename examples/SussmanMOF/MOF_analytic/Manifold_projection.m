clear
clc


hold on 
for i = 0.1:0.1:1
	interval1 = [0,1/3,0,1/3];
	interval2 = [1/3,1/2-i/6,0,1/3];
	interval3 = [1/2-i/6,1/2-i/(12-6*i),0,(1-i/2)/2+i^2/(24-12*i)];
	interval4 = [(1-i/2)/2+i^2/(24-12*i),1/2,1/2-i/(12-6*i),1/2];
	interval5 = [0,1/3,1/3,1/2-i/6];
	interval6 = [0,(1-i/2)/2+i^2/(24-12*i),1/2-i/6,1/2-i/(12-6*i)];
	interval7 = [1/2-i/(12-6*i),1/2,(1-i/2)/2+i^2/(24-12*i),1/2];
	
	fimplicit(@(x,y) y-i*x, interval1,'-.r');
	fimplicit(@(x,y) y-i/(12-24*x)-(1-2*x)*i/4,interval2,'-.b');
	fimplicit(@(x,y) y-i/(12-24*x)-(1-2*x)*i/4,interval3,'-.m');
	fimplicit(@(x,y) 1/2-(1-i)/2/((1-x)-(1-y)*i)*(1-x)-1/3*sqrt(2*i*((2*(1-x)-2*i*(1-y)-1+i)/(2*(1-x)-2*i*(1-y)))^3), interval4,'-.k');
	fimplicit(@(x,y) y-1/i*x, interval1,'-.r');
	fimplicit(@(x,y) x-i/(12-24*y)-(1-2*y)*i/4,interval5,'-.b');
	fimplicit(@(x,y) x-i/(12-24*y)-(1-2*y)*i/4,interval6,'-.m');
	fimplicit(@(x,y) (1-i)/2/((1-y)-i*(1-x))*(1-y)-1/2+1/3*sqrt(i*2*(1-(1-i)/2/((1-y)-i*(1-x)))^3), interval7,'-.k');

	interval1(1:2) = 1-interval1(1:2);k = interval1(1);interval1(1)=interval1(2);interval1(2)=k;
	interval2(1:2) = 1-interval2(1:2);k = interval2(1);interval2(1)=interval2(2);interval2(2)=k;
	interval3(1:2) = 1-interval3(1:2);k = interval3(1);interval3(1)=interval3(2);interval3(2)=k;
	interval4(1:2) = 1-interval4(1:2);k = interval4(1);interval4(1)=interval4(2);interval4(2)=k;
	interval5(1:2) = 1-interval5(1:2);k = interval5(1);interval5(1)=interval5(2);interval5(2)=k;
	interval6(1:2) = 1-interval6(1:2);k = interval6(1);interval6(1)=interval6(2);interval6(2)=k;
	interval7(1:2) = 1-interval7(1:2);k = interval7(1);interval7(1)=interval7(2);interval7(2)=k;
	fimplicit(@(x,y) y-i*(1-x), interval1,'-.r');
	fimplicit(@(x,y) y-i/(12-24*(1-x))-(1-2*(1-x))*i/4,interval2,'-.b');
	fimplicit(@(x,y) y-i/(12-24*(1-x))-(1-2*(1-x))*i/4,interval3,'-.m');
	fimplicit(@(x,y) 1/2-(1-i)/2/((1-(1-x))-(1-y)*i)*(1-(1-x))-1/3*sqrt(2*i*((2*(1-(1-x))-2*i*(1-y)-1+i)/(2*(1-(1-x))-2*i*(1-y)))^3), interval4,'-.k');
	fimplicit(@(x,y) y-1/i*(1-x), interval1,'-.r');
	fimplicit(@(x,y) (1-x)-i/(12-24*y)-(1-2*y)*i/4,interval5,'-.b');
	fimplicit(@(x,y) (1-x)-i/(12-24*y)-(1-2*y)*i/4,interval6,'-.m');
	fimplicit(@(x,y) (1-i)/2/((1-y)-i*(1-(1-x)))*(1-y)-1/2+1/3*sqrt(i*2*(1-(1-i)/2/((1-y)-i*(1-(1-x))))^3), interval7,'-.k');

	interval1(3:4) = 1-interval1(3:4);k = interval1(3);interval1(3)=interval1(4);interval1(4)=k;
	interval2(3:4) = 1-interval2(3:4);k = interval2(3);interval2(3)=interval2(4);interval2(4)=k;
	interval3(3:4) = 1-interval3(3:4);k = interval3(3);interval3(3)=interval3(4);interval3(4)=k;
	interval4(3:4) = 1-interval4(3:4);k = interval4(3);interval4(3)=interval4(4);interval4(4)=k;
	interval5(3:4) = 1-interval5(3:4);k = interval5(3);interval5(3)=interval5(4);interval5(4)=k;
	interval6(3:4) = 1-interval6(3:4);k = interval6(3);interval6(3)=interval6(4);interval6(4)=k;
	interval7(3:4) = 1-interval7(3:4);k = interval7(3);interval7(3)=interval7(4);interval7(4)=k;

	fimplicit(@(x,y) (1-y)-i*(1-x), interval1,'-.r');
	fimplicit(@(x,y) (1-y)-i/(12-24*(1-x))-(1-2*(1-x))*i/4,interval2,'-.b');
	fimplicit(@(x,y) (1-y)-i/(12-24*(1-x))-(1-2*(1-x))*i/4,interval3,'-.m');
	fimplicit(@(x,y) 1/2-(1-i)/2/((1-(1-x))-(1-(1-y))*i)*(1-(1-x))-1/3*sqrt(2*i*((2*(1-(1-x))-2*i*(1-(1-y))-1+i)/(2*(1-(1-x))-2*i*(1-(1-y))))^3), interval4,'-.k');
	fimplicit(@(x,y) (1-y)-1/i*(1-x), interval1,'-.r');
	fimplicit(@(x,y) (1-x)-i/(12-24*(1-y))-(1-2*(1-y))*i/4,interval5,'-.b');
	fimplicit(@(x,y) (1-x)-i/(12-24*(1-y))-(1-2*(1-y))*i/4,interval6,'-.m');
	fimplicit(@(x,y) (1-i)/2/((1-(1-y))-i*(1-(1-x)))*(1-(1-y))-1/2+1/3*sqrt(i*2*(1-(1-i)/2/((1-(1-y))-i*(1-(1-x))))^3), interval7,'-.k');

	interval1(1:2) = 1-interval1(1:2);k = interval1(1);interval1(1)=interval1(2);interval1(2)=k;
	interval2(1:2) = 1-interval2(1:2);k = interval2(1);interval2(1)=interval2(2);interval2(2)=k;
	interval3(1:2) = 1-interval3(1:2);k = interval3(1);interval3(1)=interval3(2);interval3(2)=k;
	interval4(1:2) = 1-interval4(1:2);k = interval4(1);interval4(1)=interval4(2);interval4(2)=k;
	interval5(1:2) = 1-interval5(1:2);k = interval5(1);interval5(1)=interval5(2);interval5(2)=k;
	interval6(1:2) = 1-interval6(1:2);k = interval6(1);interval6(1)=interval6(2);interval6(2)=k;
	interval7(1:2) = 1-interval7(1:2);k = interval7(1);interval7(1)=interval7(2);interval7(2)=k;

	fimplicit(@(x,y) (1-y)-i*x, interval1,'-.r');
	fimplicit(@(x,y) (1-y)-i/(12-24*x)-(1-2*x)*i/4,interval2,'-.b');
	fimplicit(@(x,y) (1-y)-i/(12-24*x)-(1-2*x)*i/4,interval3,'-.m');
	fimplicit(@(x,y) 1/2-(1-i)/2/((1-x)-(1-(1-y))*i)*(1-x)-1/3*sqrt(2*i*((2*(1-x)-2*i*(1-(1-y))-1+i)/(2*(1-x)-2*i*(1-(1-y))))^3), interval4,'-.k');
	fimplicit(@(x,y) (1-y)-1/i*x, interval1,'-.r');
	fimplicit(@(x,y) x-i/(12-24*(1-y))-(1-2*(1-y))*i/4,interval5,'-.b');
	fimplicit(@(x,y) x-i/(12-24*(1-y))-(1-2*(1-y))*i/4,interval6,'-.m');
	fimplicit(@(x,y) (1-i)/2/((1-(1-y))-i*(1-x))*(1-(1-y))-1/2+1/3*sqrt(i*2*(1-(1-i)/2/((1-(1-y))-i*(1-x)))^3), interval7,'-.k');

end

interval8 = [1/3,2/3,1/3,2/3];
fimplicit(@(x,y) x-y,interval8,'-.k')
interval8 = [1/3,2/3,1/3,2/3];
fimplicit(@(x,y) x+y-1,interval8,'-.k')

fimplicit(@(x,y) y-0.5*x, interval1,'r')

fimplicit(@(x,y) 1/2-x-0.5/(24*y-0.24),interval2)

% (1-x)

for V=0.05:0.05:0.5
	interval1 = [0,1/3,0,1/3];
	interval2 = [1/3,1/2,0,1/3];
	interval3 = [0,1/3,1/3,1/2];
	interval4 = [1/3,1/2,0,(1-V)/2+1/24/(1-V)*(2*(1-V)-2)^2];
	interval5 = [0,(1-V)/2+1/24/(1-V)*(2*(1-V)-2)^2,1/3,1/2];
	interval6 = [(1-V)/2+1/24/(1-V)*(2*(1-V)-2)^2,1/2,(1-V)/2+1/24/(1-V)*(2*(1-V)-2)^2,1/2];
	fimplicit(@(x,y) 9/2*x*y-V, interval1,'r');
	fimplicit(@(x,y) y-V/2-6*V*(0.5-x)^2, interval2,'b');
	fimplicit(@(x,y) x-V/2-6*V*(0.5-y)^2, interval3,'b');
	fimplicit(@(x,y) y-(1-V)/2-6*(1-V)*(0.5-x)^2, interval4,'m');
	fimplicit(@(x,y) x-(1-V)/2-6*(1-V)*(0.5-y)^2, interval5,'m');
	fimplicit(@(x,y) (1-V)^2*(4.5*(1-x)*(1-y)-3)+(1-V)*(-2.25*(1-x)-2.25*(1-y)+3)+1/8+(1-V)^3, interval6,'k');

	interval1(1:2) = 1-interval1(1:2);k = interval1(1);interval1(1)=interval1(2);interval1(2)=k;
	interval2(1:2) = 1-interval2(1:2);k = interval2(1);interval2(1)=interval2(2);interval2(2)=k;
	interval3(1:2) = 1-interval3(1:2);k = interval3(1);interval3(1)=interval3(2);interval3(2)=k;
	interval4(1:2) = 1-interval4(1:2);k = interval4(1);interval4(1)=interval4(2);interval4(2)=k;
	interval5(1:2) = 1-interval5(1:2);k = interval5(1);interval5(1)=interval5(2);interval5(2)=k;
	interval6(1:2) = 1-interval6(1:2);k = interval6(1);interval6(1)=interval6(2);interval6(2)=k;

	fimplicit(@(x,y) 9/2*(1-x)*y-V, interval1,'r');
	fimplicit(@(x,y) y-V/2-6*V*(0.5-(1-x))^2, interval2,'b');
	fimplicit(@(x,y) (1-x)-V/2-6*V*(0.5-y)^2, interval3,'b');
	fimplicit(@(x,y) y-(1-V)/2-6*(1-V)*(0.5-(1-x))^2, interval4,'m');
	fimplicit(@(x,y) (1-x)-(1-V)/2-6*(1-V)*(0.5-y)^2, interval5,'m');
	fimplicit(@(x,y) (1-V)^2*(4.5*(1-(1-x))*(1-y)-3)+(1-V)*(-2.25*(1-(1-x))-2.25*(1-y)+3)+1/8+(1-V)^3, interval6,'k');

	interval1(3:4) = 1-interval1(3:4);k = interval1(3);interval1(3)=interval1(4);interval1(4)=k;
	interval2(3:4) = 1-interval2(3:4);k = interval2(3);interval2(3)=interval2(4);interval2(4)=k;
	interval3(3:4) = 1-interval3(3:4);k = interval3(3);interval3(3)=interval3(4);interval3(4)=k;
	interval4(3:4) = 1-interval4(3:4);k = interval4(3);interval4(3)=interval4(4);interval4(4)=k;
	interval5(3:4) = 1-interval5(3:4);k = interval5(3);interval5(3)=interval5(4);interval5(4)=k;
	interval6(3:4) = 1-interval6(3:4);k = interval6(3);interval6(3)=interval6(4);interval6(4)=k;

	fimplicit(@(x,y) 9/2*(1-x)*(1-y)-V, interval1,'r');
	fimplicit(@(x,y) (1-y)-V/2-6*V*(0.5-(1-x))^2, interval2,'b');
	fimplicit(@(x,y) (1-x)-V/2-6*V*(0.5-(1-y))^2, interval3,'b');
	fimplicit(@(x,y) (1-y)-(1-V)/2-6*(1-V)*(0.5-(1-x))^2, interval4,'m');
	fimplicit(@(x,y) (1-x)-(1-V)/2-6*(1-V)*(0.5-(1-y))^2, interval5,'m');
	fimplicit(@(x,y) (1-V)^2*(4.5*(1-(1-x))*(1-(1-y))-3)+(1-V)*(-2.25*(1-(1-x))-2.25*(1-(1-y))+3)+1/8+(1-V)^3, interval6,'k');

	interval1(1:2) = 1-interval1(1:2);k = interval1(1);interval1(1)=interval1(2);interval1(2)=k;
	interval2(1:2) = 1-interval2(1:2);k = interval2(1);interval2(1)=interval2(2);interval2(2)=k;
	interval3(1:2) = 1-interval3(1:2);k = interval3(1);interval3(1)=interval3(2);interval3(2)=k;
	interval4(1:2) = 1-interval4(1:2);k = interval4(1);interval4(1)=interval4(2);interval4(2)=k;
	interval5(1:2) = 1-interval5(1:2);k = interval5(1);interval5(1)=interval5(2);interval5(2)=k;
	interval6(1:2) = 1-interval6(1:2);k = interval6(1);interval6(1)=interval6(2);interval6(2)=k;

	fimplicit(@(x,y) 9/2*x*(1-y)-V, interval1,'r');
	fimplicit(@(x,y) (1-y)-V/2-6*V*(0.5-x)^2, interval2,'b');
	fimplicit(@(x,y) x-V/2-6*V*(0.5-(1-y))^2, interval3,'b');
	fimplicit(@(x,y) (1-y)-(1-V)/2-6*(1-V)*(0.5-x)^2, interval4,'m');
	fimplicit(@(x,y) x-(1-V)/2-6*(1-V)*(0.5-(1-y))^2, interval5,'m');
	fimplicit(@(x,y) (1-V)^2*(4.5*(1-x)*(1-(1-y))-3)+(1-V)*(-2.25*(1-x)-2.25*(1-(1-y))+3)+1/8+(1-V)^3, interval6,'k');

end
