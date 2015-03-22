%Perform summation given in problem to calculate temperature 
%at a given co-ordinate

function[val] = summation(x, y, n, q)
	val = 0;
	for i in 1:n
		val = val + 2*q*((-1)^(n+1)+1)*sin(n*pi*x)*sinh(n*pi*y)/(n^2*pi^2*cosh(n*pi));
	end
end

function[] = part_b()
	%dimensions of the plate
	x = 20;
	y = 20;
	%number of steps to do in the summation
	n = 10;
	%heat flow constant
	q = .05;

	%initialize temperature matrix
	T = zeros(x,y);
	%fill in the values for the temperature matrix
	for i in 1:x
		for j in 1:y
			T(i,j) = summation(i, j, n, q)
		end
	end

	%display T
	T
end