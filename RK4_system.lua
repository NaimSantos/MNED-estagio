function f1(x,y,z)
	return -0.5*y;
end
function f2(x,y,z)
	return 4 - 0.3*z - 0.1*y;
end

local xi,xf,h = 0, 2, 0.1
local Ys = {4, 6}
local n = math.floor((xf - xi) / h)

local k = {}
for i=1,5 do
	k[i] = {}
end
for i=1,n do
	k[1][1] = f1(xi, Ys[1], Ys[2])
	k[2][1] = f2(xi, Ys[1], Ys[2])

	k[1][2] = f1(xi + 0.5*h, Ys[1] + 0.5*k[1][1]*h, Ys[2] + 0.5*k[2][1]*h)
	k[2][2] = f2(xi + 0.5*h, Ys[1] + 0.5*k[1][1]*h, Ys[2] + 0.5*k[2][1]*h)

	k[1][3] = f1(xi + 0.5*h, Ys[1] + 0.5*k[1][2]*h, Ys[2] + 0.5*k[2][2]*h)
	k[2][3] = f2(xi + 0.5*h, Ys[1] + 0.5*k[1][2]*h, Ys[2] + 0.5*k[2][2]*h)

	k[1][4] = f1(xi + h, Ys[1] + k[1][2]*h, Ys[2] + k[2][2]*h)
	k[2][4] = f2(xi + h, Ys[1] + k[1][2]*h, Ys[2] + k[2][2]*h)

	k[1][5] = (1/6)*(k[1][1] + 2*k[1][2] + 2*k[1][3] + k[1][4])*h
	k[2][5] = (1/6)*(k[2][1] + 2*k[2][2] + 2*k[2][3] + k[2][4])*h

	Ys[1] = Ys[1] + k[1][5]
	Ys[2] = Ys[2] + k[2][5]
	xi = xi + h

	print(i, xi, Ys[1], Ys[2])

end