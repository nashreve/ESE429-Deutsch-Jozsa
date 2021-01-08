
% ESE 429 
% Final Project
format compact

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                         DEUTSCH JOZSA ALGORITHM                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Modify this n and f accordingly. "n" is the number of qubits, and f is
% either a constant or balanced function

n = 2
f = @(x) rem(x,2) % This finds the remainder after dividing by 2, which is
                  % either 0 or 1. Therefore, this function is balanced. A
                  % constant function could either be "0" or "1" in place
                  % of "rem(x,2)". 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                   STEP 1: Prepping the initial state                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start by preparing p0, which is a single vector that consists of two
% different vectors: The first is an n-qubit vector initialized to |0>, and
% the second is a 1-qubit vector initialized to |1>. Thus, we use 2^n for
% the first vector, but since p0 is a combination of both, we use 2^(n+1)
% to include the second vector. Set the last bit equal to zero to load |1>
% into the second vector.

p0 = zeros(2^(n+1), 1);
p0(2) = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                   STEP 2: Create the Hadamard Matrix                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create the Hn matrix, which will allow us to apply the Hadamard Gate to
% each qubit. We use the kron function to return the Kroenecker Tensor
% Product of Hn and the Hadamard Gate. (Show each iteration as we build the
% Hn matrix). To have n Hadamard Gates in parallel is taking the tensor
% product of n gates. 

Hadamard = 1 / sqrt(2) * [1 1; 1 -1];
Hn = 1;
for i = 1:n+1
    Hn = kron(Hn, Hadamard);
    %i
end
%Hn


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                   STEP 3: Create the Oracle Function                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create U_f, which is our oracle. We have decided to construct U_f as a
% unitary matrix, so that it can perform our function on each of our
% qubits. This is because U_f is a block diagonal matrix. Each block in the
% block diagonal matrix will operate on a different qubit. There are n+1
% blocks for n+1 qubits.

% A balanced algorithm will have phase kickback add a negative phase to
% exactly half of the states. A constant algorithm will not do that. 

% To construct this U_f matrix, we use an exchange matrix "A", which is an
% antidiagonal matrix. We then build a diagonal matrix out of these "A"
% matrices, but raise "A" to the value of f(i) for i between 1 and n+1.
% Since f(i) can only be 1 or 0, this means that A will either be:

% A = [0 1; 1 0]^1 = [0 1; 1 0] (This is the exchange matrix)
% A = [0 1; 1 0]^0 = [1 0; 0 1] (This is the identity matrix)

% Any qubit multiplied by the exchange matrix will be flipped, and any
% qubit multiplied by the identity matrix will remain the same. Thus, a
% balanced function will yield exactly half of the "A" matrices in U_f to
% be exchange, and half to be identity. A constant funciton will only raise
% "A" to either 0 or 1, and thus all of the "A" matrices in U_f will be
% either exchange OR identity. This means that for a constant function the
% phases will either all remain the same (when multiplied by identity), or
% they will all flip (when multiplied by exchange- but this is essentially
% the same as not flipping in the first place since everyone is flipping).

A = [0 1; 1 0];
U_f = blkdiag(A^f(0), A^f(1));
for i = 2:2^n-1
    U_f = blkdiag(U_f, A^f(i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                   STEP 4: Running the algorithm                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Now, we simulate the algorithm at each of the different stages. p1 is
% simply taking the Hadamard of each of the initial qubits. p2 is running
% the qubits through the oracle. p3 is taking the Hadamard of the oracle's
% output in order to bring the states back down from superpositon. 

P0 = p0.'
p1 = Hn * p0; P1 = p1.'
p2 = U_f * p1; P2 = p2.'
p3 = Hn * p2; P3 = p3.'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                 STEP 5: Evaluating whether CONST or BAL                 %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Now, measure all elements of p3 except for the last two. This is because
% the last two elemnents represent the first n-qubit vector that was
% intialized to |0> which is what we care about. 

% If they're all zero, then f is CONSTANT. If any of them are non-zero,
% then the result is BALANCED. Start by assuming "constant", then switch to
% "balanced" if we find a non-zero after iterating through p3. 

result = ' CONSTANT';
for i = 3:(length(p3))
    rounded = round(p3(i),4);
    if (rounded ~= 0) || (rounded ~= (-0))    
        result = 'BALANCED';
    end
end
result

figure;
imagesc(U_f);
s = strcat(func2str(f), ' is ' , result);
title(s);
colorbar;


