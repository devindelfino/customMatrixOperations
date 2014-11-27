%Filename: Matrix.m
%Description: Implements a custom Matrix class with standard matrix operations (see README.md for more details)

classdef Matrix
	properties(SetAccess = private)
		r = 2;
		c = 2;
		data = zeros(2);
	end %end properties

	methods
		% CONSTRUCTOR ------------------------------------------------------------------------------------------
		function M = Matrix(mat)
		    % Description: default constructor
		    % Parameters: mat - a standard MATLAB matrix
		    % Returns: an instance of the custom Matrix class
			if(nargin == 0)
				M.r = 2;
				M.c = 2;
				data = zeros(2);
			end
			temp_size = size(mat);
			M.r = temp_size(1);
			M.c = temp_size(2);
			M.data = mat;
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

		function elmt = get_elements(A, i, j)
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
						diff_mat(itR, itC) = A.get_element(itR, itC) - B.get_element(itR, itC);
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
						size(A.get_elements(itR,:))
						size(B.get_elements(:,itC))
						prod_mat(itR, itC) = sum(A.get_elements(itR,:) .* transpose(B.get_elements(:,itC)));
					end
				end

				M = Matrix(prod_mat);

			else % the number of columns in A and the number of columns in B don't agree.
				% raise error
				error('Inner Matrix dimensions must agree (the number of columns in the first matrix must equal the number of rows in the second matrix)')
			end
		end

		% MATRIX INVERSION ------------------------------------------------------------------------------------------
		function M = invert(A) 
			% Description: Calculates the inverse matrix of Matrix A
			% Parameters: A - a custom Matrix object
			% Returns: The inverse matrix or [0] if not invertible 

			
		end

		% DETERMINANT ------------------------------------------------------------------------------------------
		function determinant = det(A) 
			% Description: Calculates the determinant of matrix A
			% Parameters: A - a custom Matrix object
			% Returns: The determinant of Matrix A

			
		end
	end % end methods
end % end classdef
