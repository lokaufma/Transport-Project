%Perform summation given in problem to calculate temperature 
%at a given co-ordinate

function[] = part_b()
	%dimensions of the plate
	w = 20;
	l = 20;
	%number of steps to do in the summation
	n = 10;
	%heat flow constant
	q = .05;

	%initialize temperature matrix
	T = zeros(w,l);
	%fill in the values for the temperature matrix
	for i = 1:w
		for j = 1:l
			T(i,j) = summation(i, j, n, q, w, l);
		end
	end

	%display T
	imagesc(T);
end

function[val] = summation(x, y, n, q, w, l)
	val = 0;
	for i = 1:n
		val = val + 2*q*((-1)^(i+1)+1)*sin(i*pi*x/w)*sinh(i*pi*y/l)/(i^2*pi^2*cosh(i*pi));
	end
end