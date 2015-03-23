% Perform summation given in problem to calculate temperature 
% at a given co-ordinate

function[] = part_b()
	% Dimensions of the plate
	w = 25;
	l = 25;
	% Number of steps to do in the summation
	n = 100;
	% Heat flow constant
	q = .05;

	% Initialize temperature matrix
	T = zeros(w,l);
	% Fill in the values for the temperature matrix
	for i = 1:w
		for j = 1:l
			T(j, i) = summation(i, j, n, q, w, l);
		end
	end

	% Display T
	imagesc(T);
    set(gca, 'YDir', 'normal');
    colorbar;
end

function[val] = summation(x, y, n, q, w, l)
	val = 0;
	for i = 1:n
		val = val + 2*q*((-1)^(i+1)+1)*sin(i*pi*x/w)*sinh(i*pi*y/l)/(i^2*pi^2*cosh(i*pi));
	end
end