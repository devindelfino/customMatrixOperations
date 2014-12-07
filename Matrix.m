%Filename: Matrix.m
%Description: Implements a custom Matrix class with standard matrix operations (see README.md for more details)

classdef Matrix
	properties(SetAccess = private)
		r = 2;
		c = 2;
		data = zeros(2);
	end %end properties (private)

	properties
		determinant = NaN;
		invertible = 0;
		dominant_eigen_value = NaN;
	end

	methods
		% CONSTRUCTOR ------------------------------------------------------------------------------------------
		function M = Matrix(mat)
		    % Description: default constructor
		    % Parameters: mat - a standard MATLAB matrix
		    % Returns: an instance of the custom Matrix class
			if(nargin == 0)
				M.r = 2;
				M.c = 2;
				M.data = zeros(2);
				M.determinant = NaN;
				M.invertible = 0;
			end
			temp_size = size(mat);
			M.r = temp_size(1);
			M.c = temp_size(2);
			M.data = mat;

			if(M.r == M.c)
				if(det(mat) == 0)
					M.invertible = 0;
					M.determinant = 0;
				else
					M.invertible = 1;
				end
			else
				M.invertible = 0;
				M.determinant = NaN;
			end
			M.dominant_eigen_value = 0;
		end

		% DISPLAY ------------------------------------------------------------------------------------------
		function display(A)
			% Description: Displays the contents of the Matrix
			% Parameters: None
			% Returns: A standard MATLAB matrix
			A.data
		end

		% GET DIMENSIONS ------------------------------------------------------------------------------------------
		function dim = dims(A)
			% Description: Prints the matrix dimensions
			% Parameters: None
			% Returns: A 1x2 matrix where the first element is the number of rows and the second element is the number of columns 
			% display([num2str(M.r) ' rows by ' num2str(M.c) ' columns']);
			dim = [A.r A.c];
		end

		% GET ELEMENTS ------------------------------------------------------------------------------------------
		function elmt = get_elements(A, i, j)
			% Description: Gets an element(s) of a Matrix
			% Parameters: i - integer or array of integers indicating the rows (: for all)
			%             j - integer or array of integers indicating the columns (: for all)
			% Returns: An element or list of elements from the Matrix A 
			elmt = A.data(i, j);
		end

		% MATRIX TRANSPOSITION  ------------------------------------------------------------------------------------------
		function M = transpose(A) % equivalent to (A.')
			% Description: Transposes the Matrix
			% Parameters: None
			% Returns: The transposition Matrix of Matrix A
			temp_size = A.dims();
			trans_mat = zeros(temp_size(2), temp_size(1));

			for(itR = 1:temp_size(2))
				trans_mat(itR,:) = A.get_elements(:,itR);
			end

			M = Matrix(trans_mat);
		end

		% MATRIX EQUALITY (custom, MATLAB) ------------------------------------------------------------------------------------------

		function tf = eq(A,b)
			% Description: Tests for equality between A and B to a factor of 0.0001
			% Parameters: A - a custom Matrix object
			%             B - a matlab Matrix object
			% Returns: boolean stating whether or not A == B
			
			ERR = 0.0001;
			adims = A.dims();
			bdims = size(b);

			a = A.get_elements(:,:);

			if(adims == bdims)
				if(all(abs(a-b) < ERR))
					tf = 1;
				else
					tf = 0;
				end
			else
				tf = 0;
			end
		end

		% MATRIX ADDITION ------------------------------------------------------------------------------------------

		function M = plus(A, B) % equivalent to (A + B)
			% Description: Adds Matrix A and Matrix B (A and B must have the same dimensions)
			% Parameters: A - a custom Matrix object
			%             B - a custom Matrix object
			% Returns: The matrix containing the sum of matrices A and B

			if(A.dims() == B.dims())
				temp_size = A.dims();
				sum_mat = zeros(temp_size(1), temp_size(2));

				for(itR = 1:temp_size(1))
					for(itC = 1:temp_size(2))
						sum_mat(itR, itC) = A.get_elements(itR, itC) + B.get_elements(itR, itC);
					end
				end

				M = Matrix(sum_mat);

			else % matrix dimensions don't agree
				% raise error
				error('Matrix dimensions must agree.')
			end
		end

		% MATRIX SUBTRACTION ------------------------------------------------------------------------------------------
		function M = minus(A, B) % equivalent to (A - B)
			% Description: Subtracts Matrix B from Matrix A (A and B must have the same dimensions)
			% Parameters: A - a custom Matrix object
			%             B - a custom Matrix object
			% Returns: The matrix containing the sum of matrices A and B 

			if(A.dims() == B.dims())
				temp_size = A.dims();
				diff_mat = zeros(temp_size(1), temp_size(2));

				for(itR = 1:temp_size(1))
					for(itC = 1:temp_size(2))
						diff_mat(itR, itC) = A.get_elements(itR, itC) - B.get_elements(itR, itC);
					end
				end

				M = Matrix(diff_mat);

			else % matrix dimensions don't agree
				% raise error
				error('Matrix dimensions must agree.')
			end
		end

		% MATRIX MULTIPLICATION ------------------------------------------------------------------------------------------
		function M = mtimes(A, B) % equivalent to (A * B)
			% Description: Matrix Multiplies Matrix A and Matrix A (A's columns must equal B's rows)
			% Parameters: A - a custom Matrix object
			%             B - a custom Matrix object
			% Returns: The matrix containing the product of matrices A and B 

			A_size = A.dims();
			B_size = B.dims();

			if(A_size(2) == B_size(1))
				prod_mat = zeros(A_size(1), B_size(2));
				for(itR = 1:A_size(1))
					for(itC = 1:B_size(2))
						% size(A.get_elements(itR,:))
						% size(B.get_elements(:,itC))
						prod_mat(itR, itC) = sum(A.get_elements(itR,:) .* transpose(B.get_elements(:,itC)));
					end
				end

				M = Matrix(prod_mat);

			else % the number of columns in A and the number of columns in B don't agree.
				% raise error
				error('Inner Matrix dimensions must agree (the number of columns in the first matrix must equal the number of rows in the second matrix)')
			end
		end

		% LU FACTORIZATION WITH PARTIAL PIVOTING ------------------------------------------------------------------------------------------
		% Used Algorithm from http://www.math.iit.edu/~fass/477577_Chapter_7.pdf, page 66
		function [L U P]= LU_factor(A)
			% Description: Calculates the Upper and Lower triangular matrices of square Matrix A along with the permutation Matrix P for stability
			% Parameters: A - a custom Matrix object
			% Returns: The Upper triangular Matrix (U), the Lower triangular Matrix (L), and the Permutation Matrix (P)
			dimensions = A.dims();
			m = dimensions(1); n = dimensions(2);

			if(m == n) % check for a square matrix
				perm_parity = 0; %used to predetermine determinant of permuation matrix
				Upper = A.get_elements(:,:); Lower = eye(m); Perm = eye(m);
				for(k = 1:(m-1))

					for(i = k:m)
						swap = Upper(k, k:m);
						Upper(k, k:m) = Upper(i, k:m);
						Upper(i, k:m) = swap;

						swap = Lower(k, 1:k-1);
						Lower(k, 1:k-1) = Lower(i, 1:k-1);
						Lower(i, 1:k-1) = swap;

						swap = Perm(k, :);
						Perm(k, :) = Perm(i, :);
						Perm(i, :) = swap;
						if(i ~= k)
							perm_parity = perm_parity + 1;
						end
					end

					for(j = (k+1):m)
						if(Upper(k, k) == 0)
							Lower(j, k) = 0;
						else
							Lower(j, k) = Upper(j, k) / Upper(k, k);
						end
						Upper(j, k:m) = Upper(j, k:m) - Lower(j, k) * Upper(k, k:m);
					end
				end


				U = Matrix(Upper); L = Matrix(Lower); P = Matrix(Perm);

				if(mod(perm_parity, 2) == 0)
					P.determinant = 1;
				else
					P.determinant = -1;
				end

			else
				error('Starting Matrix must be square.')
			end
		end

		% DETERMINANT ------------------------------------------------------------------------------------------
		function determinant = det(A) 
			% Description: Calculates the determinant of matrix A (if already calculated, returns determinant)
			% Parameters: A - a custom Matrix object
			% Returns: The determinant of Matrix A
			dimensions = A.dims();
			m = dimensions(1); n = dimensions(2);

			if(m == n) % check for square matrix
				
				if(isnan(A.determinant))
					[L U P] = A.LU_factor();

					detP = P.determinant;
					detL = 1;
					detU = 1;

					for(i = 1:m) % calculating the product of the diagonal entries of Lower and Upper triangular matrices of A
						detL = detL * L.get_elements(i,i);
						detU = detU * U.get_elements(i,i);
					end

					determinant = detP * detL * detU; % A = P * L * U => det(A) = det(P) * det(L) * det(U)
					A.determinant = determinant;
				else
					determinant = A.determinant;
				end

			else
				error('Starting Matrix must be square.')
			end
		end

		% NORMALIZE ------------------------------------------------------------------------------------------
		function M = normalize(A)
			% Description: Normalizes the vector A
			% Parameters: A -  Vector A
			% Returns: The normalized vector M
			dimensions = A.dims();
			m = dimensions(1); n = dimensions(2);

			if(m == 1 | n == 1)
				temp = A.get_elements(:,:);
				len = sqrt(sum(temp.^2));
				M = Matrix(temp./len);
			else% matrix dimensions don't agree
				% raise error
				error('Only vectors can be normalized (one column or one row).')
			end
		end

		% MAX ------------------------------------------------------------------------------------------
		function maximum = max(A) 
			% Description: Finds the maximum element in a Vector
			% Parameters: A - a custom Matrix object
			% Returns: The maximum element
			
			maximum = max(A.data);
		end

		% SCALAR DIVISION ------------------------------------------------------------------------------------------
		function M = rdivide(A, scalar) 
			% Description: Divides each element in the Matrix by a scalar (equivalent to A ./ scalar)
			% Parameters: A - a custom Matrix object
			%             scalar - a scalar
			% Returns: The quotient matrix
			
			if(scalar == 0)
				error('Cannot divide by zero');
			else
				M = Matrix(A.data ./ scalar);
			end
		end

		% 1x1 MATRIX DIVISION ------------------------------------------------------------------------------------------
		function quotient = mrdivide(A, B) 
			dimensions = A.dims();
			mA = dimensions(1); nA = dimensions(2);
			dimensions = B.dims();
			mB = dimensions(1); nB = dimensions(2);

			if(mA == 1 & nA == 1 & mB == 1 & nB == 1)
				quotient = A.data ./ B.data;
			else
				error('Can only divide single values');
			end
		end

		% EIGENVALUES (POWER TECHNIQUE) ------------------------------------------------------------------------------------------
		% Source: http://www.robots.ox.ac.uk/~sjrob/Teaching/EngComp/linAlg34.pdf, page 81
		function [eig_vector eig_value] = eig_dominant(A, iterations) 
			% Description: Calculates the dominant eigenvalue of the matrix using the Power Method Iterative Technique
			% Parameters: A - a custom Matrix object
			%             iterations - an integer representing the number of iterations desired
			% Returns: The dominant eigenvalue	
			dimensions = A.dims();
			m = dimensions(1); n = dimensions(2);

			if(m == n)
				current_vector = Matrix(zeros(n, 1) + 1);
				next_vector = current_vector;

				for(it = 1:iterations)
					next_it = A * current_vector;
					next_vector = next_it.normalize();
					current_vector = next_vector;
				end

				eig_vector = current_vector ./ max(current_vector);
				eig_value = (eig_vector.' * (A * eig_vector )) / (eig_vector.' * eig_vector);
			else
				error('Starting Matrix must be square.')
			end
		end

		% EIGENVALUES (INVERSE TECHNIQUE) ------------------------------------------------------------------------------------------
		% Source: http://mathreview.uwaterloo.ca/archive/voli/1/panju.pdf
		function [eig_vectors eig_values] = eig_all(A, iterations) 
			% Description: Calculates the eigenvectors/eigenvalues of the matrix A using the Arnoldi Iteration Method
			% Parameters: A - a custom Matrix object
			%             iterations - an integer representing the number of iterations desired
			% Returns: All of the eigenvectors and eigenvalues of the matrix A	
			
			dimensions = A.dims();
			m = dimensions(1); n = dimensions(2);

			if(m == n)
				shift = Matrix((est*eye(m)));
				Ainv = A - shift;
				current_vector = Matrix(zeros(n, 1) + 1);
				next_vector = current_vector;

				for(it = 1:iterations)

					next_it = ge_linsolve(Ainv, current_vector);
					next_vector = next_it.normalize();
					current_vector = next_vector;
				end

				eig_vector = current_vector ./ max(current_vector);
				eig_value = (eig_vector.' * (A * eig_vector )) / (eig_vector.' * eig_vector);
			else
				error('Starting Matrix must be square.')
			end
		end

		% SOLVE LINEAR SYSTEM using LU FACTORIZATION AND GAUSSIAN ELIMINATION ------------------------------------------------------------------------------------------
		% Source: http://www.cs.mtu.edu/~shene/COURSES/cs3621/NOTES/INT-APP/CURVE-linear-system.html
		function X = ge_linsolve(A, B) 
			% Description: Solve the linear system Ax = b for a square matrix A using LU Factorization and Gaussian Elimination (forward/backward substitution)
			% Parameters: A - an m x m custom Matrix object
			%             B - an m x k custom Matrix object
			% Returns: The linear solution as an m x k custom Matrix object
			dimensions = A.dims();
			m = dimensions(1); n = dimensions(2);

			if(m == n)
				[L U P] = A.LU_factor();

				% P * A = L * U =====> A = invP * L * U
				
				% A * X = B now becomes (L * U) * X = P * B

				% 1) Since (L * U) * X = (P * B), => L * (U * X) = (P * B) => L * Y = (P * B) where Y = U * X (SOLVE FOR Y using forward substitution)
				% 2) Since Y = U * X (SOLVE FOR X using backward substitution)

				% 1) solve for Y, L * Y = P*B
				% PB = B;
				PB = P * B;
				PBdims = PB.dims();
				Y = zeros(n,PBdims(2));

				% foward substitution
				for(c = 1:PBdims(2))	% for each column in Y
					Y(1,c) = PB.get_elements(1,c) ./ L.get_elements(1,1);

					for(r = 2:n) % elements within the column
						sum = 0;
						for(k = 1:r-1)
							sum = sum + (L.get_elements(r,k) .* Y(k,c));
						end
						Y(r,c) = (PB.get_elements(r,c) - sum) ./ L.get_elements(r,r);
					end
				end

				% 2) solve for X, Y = U * X
				
				X = zeros(n,PBdims(2));

				% backward substitution
				for(c = 1:PBdims(2))	% for each column in X
					X(n,c) = Y(n,c) ./ U.get_elements(n,n);

					for(r = n-1:-1:1) % elements within the column
						sum = 0;
						for(k = r+1:n)
							sum = sum + (U.get_elements(r,k) .* X(k,c));
						end
						X(r,c) = (Y(r,c) - sum) ./ U.get_elements(r,r);
					end
				end

				X = Matrix(X);
			else
				error('Starting Matrix must be square.')
			end
		end

		% SOLVE LINEAR SYSTEM using Gauss-Seidel Method ------------------------------------------------------------------------------------------
		function X = gs_linsolve(A, B) 
			% Description: Solve the linear system Ax = b for any matrix A using Gauss-Seidel Iterative Method
			% Parameters: A - an m x n custom Matrix object
			%             B - an m x k custom Matrix object
			% Returns: The linear solution as an n x k custom Matrix object
			dimensions = A.dims();
			m = dimensions(1); n = dimensions(2);

			
		end

		% MATRIX INVERSION ------------------------------------------------------------------------------------------
		function M = inv(A) 
			% Description: Calculates the inverse matrix of Matrix A
			% Parameters: A - a custom Matrix object
			% Returns: The inverse matrix or a matrix of inf values if not-invertible
			dimensions = A.dims();
			m = dimensions(1); n = dimensions(2);

			if(m == n) % check for square matrix
				if(A.invertible)
					% calculate inverse
					M = ge_linsolve(A, Matrix(eye(m)));
				else
					display('The matrix is singular because it has a determinant of 0 (not invertible).')
					M = zeros(m) + inf;
				end
			else
				error('Starting Matrix must be square.')
			end
		end

	end % end methods
end % end classdef
