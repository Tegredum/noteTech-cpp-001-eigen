% partialPivLUInverse_test.m
% encoding: utf-8
% author: Fitten Code Agent
% matlab version: R2024b
% Test script for partialPivLUInverse function

% Test 1: Random matrix
n = 5;
A = randn(n);
testInverse(A, 'Random matrix');

% Test 2: Diagonal matrix
A = diag([1,2,3,4,5]);
testInverse(A, 'Diagonal matrix');

% Test 3: Symmetric positive definite matrix
A = randn(n);
A = A'*A + eye(n); % make it positive definite
testInverse(A, 'Symmetric positive definite');

% Test 4: Small matrix
A = [1 2; 3 4];
testInverse(A, '2x2 matrix');

% Test 5: Hilbert matrix (ill-conditioned)
n = 4;
A = hilb(n);
testInverse(A, 'Hilbert matrix');

% Helper function to run tests
function testInverse(A, testName)
	fprintf('\n=== Testing: %s ===\n', testName);
	
	% Compute inverse using our function
	tic;
	invA = partialPivLUInverse(A);
	elapsed = toc;
	fprintf('Our function time: %.4f sec\n', elapsed);
	
	% Compute inverse using MATLAB built-in
	tic;
	invA_matlab = inv(A);
	elapsed_matlab = toc;
	fprintf('MATLAB inv time: %.4f sec\n', elapsed_matlab);
	
	% Check accuracy
	err = norm(A*invA - eye(size(A)), 'fro');
	fprintf('Error norm (A*invA-I): %.2e\n', err);
	
	% Compare with MATLAB result
	diff = norm(invA - invA_matlab, 'fro');
	fprintf('Difference from MATLAB inv: %.2e\n', diff);
	
	% Display some results
	fprintf('First element: ours=%.6f, MATLAB=%.6f\n', invA(1,1), invA_matlab(1,1));
end
