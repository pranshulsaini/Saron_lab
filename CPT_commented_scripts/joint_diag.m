function [ V , D ] =  joint_diag(A, jthresh)   % A is of dimen 88 x (88*n_points)
% function [ V , D ] =  joint_diag(A,jthresh,max_iter)
% Joint approximate diagonalization
% 
% Joint approximate of n (complex) matrices of size m*m stored in the
% m*mn matrix A by minimization of a joint diagonality criterion
%
% Usage:  [ V , D ] =  joint_diag(A,jthresh)
%
% Input :
% * the m*nm matrix A is the concatenation of n matrices with size m
%   by m. We denote A = [ A1 A2 .... An ]
% * threshold is an optional small number (typically = 1.0e-8 see the M-file).
%
% Output :
% * V is an m*m unitary matrix.
% * D = V'*A1*V , ... , V'*An*V has the same size as A and is a
%   collection of diagonal matrices if A1, ..., An are exactly jointly
%   unitarily diagonalizable.
%

% The algorithm finds a unitary matrix V such that the matrices
% V'*A1*V , ... , V'*An*V are as diagonal as possible, providing a
% kind of `average eigen-structure' shared by the matrices A1 ,...,An.
% If the matrices A1,...,An do have an exact common eigen-structure ie
% a common orthonormal set eigenvectors, then the algorithm finds it.
% The eigenvectors THEN are the column vectors of V and D1, ...,Dn are
% diagonal matrices.
% 
% The algorithm implements a properly extended Jacobi algorithm.  The
% algorithm stops when all the Givens rotations in a sweep have sines
% smaller than 'threshold'.
%
% In many applications, the notion of approximate joint
% diagonalization is ad hoc and very small values of threshold do not
% make sense because the diagonality criterion itself is ad hoc.
% Hence, it is often not necessary in applications to push the
% accuracy of the rotation matrix V to the machine precision.
%
% PS: If a numerical analyst knows `the right way' to determine jthresh
%     in terms of 1) machine precision and 2) size of the problem,
%     I will be glad to hear about it.
% 
%
% This version of the code is for complex matrices, but it also works
% with real matrices.  However, simpler implementations are possible
% in the real case.
%
% See more info, references and version history at the bottom of this
% m-file
%
%----------------------------------------------------------------
% Version 1.2
%
% Copyright 	: Jean-Francois Cardoso. 
% Author 	: Jean-Francois Cardoso. cardoso@sig.enst.fr
% Comments, bug reports, etc are welcome.
%----------------------------------------------------------------


    [m,nm] = size(A);

    % Better declare the variables used in the loop :
    B       = [ 1 0 0 ; 0 1 1 ; 0 -i i ] ;   % i means imaginary values (P)
    Bt      = B' ;
    
    % The assigned values to the variables below seem to be unused (P)
    Ip      = zeros(1,nm) ;
    Iq      = zeros(1,nm) ;
    g       = zeros(3,nm) ;  
    g	    = zeros(3,m);
    G       = zeros(2,2) ;    
    vcp     = zeros(3,3);
    D       = zeros(3,3);
    la      = zeros(3,1);
    K       = zeros(3,3);
    angles  = zeros(3,1);
    pair    = zeros(1,2);
    G	    = zeros(3);
    c       = 0 ;
    s       = 0 ;

    % Init
    V	= eye(m);
    encore	= 1; 
    iter_count=0;
    while encore,
      encore=0;            % I think this is not the right indent (P)
      smax=0;              % I think this is not the right indent (P)
      for p=1:m-1, Ip = p:m:nm ; % Ip (dim: n) is indexing through all the values corresponding to channel p in all future epochs (P)
        for q=p+1:m, Iq = q:m:nm ; % Iq (dim: n) is indexing through all the values corresponding to channel q in all future epochs (P)

          % Computing the Givens angles
          g       = [ A(p,Ip)-A(q,Iq)  ; A(p,Iq) ; A(q,Ip) ] ;  % Dim: 3 x n (P)     I don't know why this was done (P)
          [vcp,D] = eig(real(B*(g*g')*Bt));         % Dim of B*(g*g')*Bt: 3 x 3. real only extracts real values out of the matrix. vcp is eig_vec with dim: 3x3. D is eig_val with dim: 3 x3  (P)
          [la, K] = sort(diag(D));  % sorting of eigen values in ascending order (P)
          angles  = vcp(:,K(3));   % getting the eigen vector for the largest eigen value (P)
          if angles(1)<0 , angles= -angles ; end ;         % It only compares the real part of 1st element of largest eigen vector in the if condition (P)
          c       = sqrt(0.5+angles(1)/2);
          s       = 0.5*(angles(2)-j*angles(3))/c; % j here is for imaginary notation (P)

          smax=max(s,smax); % If s is imaginary, then it's norm will be compared with smax, which means it will be chosen as smax = 0 (P)

          if abs(s)>jthresh, %%% updates matrices A and V by a Givens rotation              abs(s) gives out the norm of s (P)
            encore          = 1 ;
            pair            = [p;q] ;  %  2x1 (P)
            G               = [ c -conj(s) ; s c ] ;  % 2x2 (P)
            V(:,pair)       = V(:,pair)*G ;  % 88x2  (P)
            A(pair,:)       = G' * A(pair,:) ; % 2 x mn (P)
            A(:,[Ip Iq])    = [ c*A(:,Ip)+s*A(:,Iq) -conj(s)*A(:,Ip)+c*A(:,Iq) ] ;

          end%% if_
        end%% q loop
      end%% p loop

      iter_count=iter_count+1;
      fprintf('Iteration %d Threshold: %e\n', iter_count,smax);

    end%% while_

    D = A ;  % m x nm (P)

    return

% Revision history
%
% Version 1.2.  Nov. 2, 1997.
%   o some Matlab tricks to have a cleaner code.
%   o Changed (angles=sign(angles(1))*angles) to (if angles(1)<0 ,
%   angles= -angles ; end ;) as kindly suggested by Iain Collings
%   <i.collings@ee.mu.OZ.AU>.  This is safer (with probability 0 in
%   the case of sample statistics)
%
% Version 1.1.  Jun. 97.
% 	Made the code available on the WEB




%----------------------------------------------------------------
% References:
%
% The 1st paper below presents the Jacobi trick.
% The second paper is a tech. report the first order perturbation
% of joint diagonalizers
%
%
%@article{SC-siam,
%  HTML	       = "ftp://sig.enst.fr/pub/jfc/Papers/siam_note.ps.gz",
%  author       = "Jean-Fran\c{c}ois Cardoso and Antoine Souloumiac",
%  journal      = "{SIAM} J. Mat. Anal. Appl.",
%  title 	= "Jacobi angles for simultaneous diagonalization",
%  pages 	= "161--164",
%  volume       = "17",
%  number       = "1",
%  month 	= jan,
%  year 	= {1996}
%  }
%
%
%
%@techreport{PertDJ,
%  author       = "Jean-Fran\c{c}ois Cardoso",
%  HTML	        = "ftp://sig.enst.fr/pub/jfc/Papers/joint_diag_pert_an.ps",
%  institution  = "T\'{e}l\'{e}com {P}aris",
%  title        = "Perturbation of joint diagonalizers. Ref\# 94D027",
%  year	        = "1994"
%}
end