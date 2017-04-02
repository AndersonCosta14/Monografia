function [t, x, y, z] = RungeKutta4th(EDO1, EDO2, EDO3, a, b, h, x1, y1, z1, par1, par2, par3, par4, par5, par6, par7, par8);
	%% 
	% EDO1: dx/dt
	% EDO2: dy/dt
	% EDO3: dz/dt
	% a: Primeiro Valor de t
	% b: Último Valor de t
	% h: Tamanho da Discretização
	% x1: Valor Inicial de x
	% y1: Valor Inicial de y
	% y1: Valor Inicial de z
	%%
	t(1) = a;
	x(1) = x1;
	y(1) = y1;
	z(1) = z1;

	n    = (b-a)/h;

	for i = 1:n
		t(i+1) = t(i) + h;
		tm     = t(i) + h/2;
		Kx1    = feval(EDO1, t(i), x(i), y(i), z(i), par1, par2, par3, par4, par5, par8);
		Ky1    = feval(EDO2, t(i), x(i), y(i), z(i), par1, par2, par3, par4, par5, par6, par8);
		Kz1    = feval(EDO3, t(i), x(i), y(i), z(i), par1, par2, par3, par4, par5, par7, par8);
		Kx2    = feval(EDO1, tm, x(i) + Kx1*h/2, y(i) + Ky1*h/2, z(i) + Kz1*h/2, par1, par2, par3, par4, par5, par8);
		Ky2    = feval(EDO2, tm, x(i) + Kx1*h/2, y(i) + Ky1*h/2, z(i) + Kz1*h/2, par1, par2, par3, par4, par5, par6, par8);
		Kz2    = feval(EDO3, tm, x(i) + Kx1*h/2, y(i) + Ky1*h/2, z(i) + Kz1*h/2, par1, par2, par3, par4, par5, par7, par8);
		Kx3    = feval(EDO1, tm, x(i) + Kx2*h/2, y(i) + Ky2*h/2, z(i) + Kz2*h/2, par1, par2, par3, par4, par5, par8);
		Ky3    = feval(EDO2, tm, x(i) + Kx2*h/2, y(i) + Ky2*h/2, z(i) + Kz2*h/2, par1, par2, par3, par4, par5, par6, par8);
		Kz3    = feval(EDO3, tm, x(i) + Kx2*h/2, y(i) + Ky2*h/2, z(i) + Kz2*h/2, par1, par2, par3, par4, par5, par7, par8);
		Kx4    = feval(EDO1, t(i+1), x(i) + Kx3*h, y(i) + Ky3*h, z(i) + Kz3*h, par1, par2, par3, par4, par5, par8);
		Ky4    = feval(EDO2, t(i+1), x(i) + Kx3*h, y(i) + Ky3*h, z(i) + Kz3*h, par1, par2, par3, par4, par5, par6, par8);
		Kz4    = feval(EDO3, t(i+1), x(i) + Kx3*h, y(i) + Ky3*h, z(i) + Kz3*h, par1, par2, par3, par4, par5, par7, par8);
		x(i+1) = x(i) + (Kx1 + 2*Kx2 + 2*Kx3 + Kx4)*h/6;
		y(i+1) = y(i) + (Ky1 + 2*Ky2 + 2*Ky3 + Ky4)*h/6;
		z(i+1) = z(i) + (Kz1 + 2*Kz2 + 2*Kz3 + Kz4)*h/6;
	end

end