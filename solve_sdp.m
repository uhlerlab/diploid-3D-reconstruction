function [G0, x0, cvx_status] = solve_sdp(D_ij_pairwise, D_ii_homolog, D_ijk_higher_order, indices_ijk_higher_order, D_neighboring, indices_neighboring, rho)

  % Solve SDP relaxation on the Gram matrix to obtain 3D diploid configuration.
  % The function takes as input:
  %  D_ij_pairwise: n x n matrix of distances derived from contact frequency data
  %  D_ii_homolog: 1 x n matrix of distances between homologous loci
  %  D_ijk_higher_order: n x n x n tensor of higher-order (3-way) distances between loci
  %  indices_ijk_higher_order: number of constraints x 3 matrix which stores (i, j, k) indices corresponding to the 3-way distances
  %  D_neighboring: vector of distances between neighboring beads
  %  indices_neighboring: vector of indices for neighboring beads (all beads except the last bead on the chromosome)
  %  rho: hyperparameter for weighing the trace term of the objective
  % The function returns:
  %  G0: Gram matrix, which is the solution to the SDP
  %  x0: 2n x 3 matrix corresponding to the 3D diploid configuration
  %  cvx_status: Status of CVX


	m_hap=length(D_ij_pairwise); % total number of homolog domains
	m=2*m_hap; % total number of domains
	dim=3; % dimension
	numc = length(indices_ijk_higher_order); % number of higher-order constraints


  % Define SDP
  cvx_solver mosek
  cvx_begin
  cvx_precision best
  variable lamVals(8, numc);
  variable lamSlack(8, numc);
  variable G0(m,m) symmetric;

  % Minimize the trace of the Gram matrix as an approximation to matrix rank
  val_min = trace(G0) * rho;


  % Express pairwise distances as constraints on the Gram matrix
  for i = 1:(m_hap-1)
      for j = (i+1):m_hap
          val_min = val_min + (G0(i,i)+G0(j,j)+G0(m_hap+i, m_hap+i)+G0(m_hap+j, m_hap+j)-G0(i,j)-G0(m_hap+i,j)-G0(m_hap+i, m_hap+j)-G0(i,m_hap+j) - 0.5 * D_ij_pairwise(i, j))^2;
      end;
  end;

  % Express homolog-homolog distances as constraints on the Gram matrix
  for i=1:m_hap
      val_min = val_min + ((G0(i,i)+G0(m_hap+i,m_hap+i)-2*G0(i,m_hap+i) - D_ii_homolog(i))^2);
  end;


  % Express distances between neighboring beads as constraints on the Gram matrix
  for index=1:length(indices_neighboring)
      i = indices_neighboring(index);
      i_next = indices_neighboring(index) + 1;
      dist_neighbor = D_neighboring(index);

      val_min = val_min + (G0(i,i) + G0(i_next, i_next)  - 2* G0(i, i_next) - dist_neighbor)^2;
      val_min = val_min + (G0(m_hap+i, m_hap+i) + G0(m_hap+i_next, m_hap+i_next)  - 2* G0(m_hap+i, m_hap+i_next) - dist_neighbor)^2;
  end


  % Express distances of order 3 as constraints on the Gram matrix
  val_min = val_min + sum(sum(lamSlack));
  val_min = val_min + norm(sum(lamVals, 1), 2);


  for i1 = 1:numc
      i = indices_ijk_higher_order(i1, 1);
      j = indices_ijk_higher_order(i1, 2);
      k = indices_ijk_higher_order(i1, 3);
      D_ijk_higher_order(i,j,k) + lamVals(1, i1) == lamSlack(1, i1) + G0(i, i) + G0(j, j) + G0(k, k) -  G0(i, j) -  G0(i, k) - G0(j, k);
      D_ijk_higher_order(i,j,k) + lamVals(2, i1) == lamSlack(2, i1) + G0(i, i) + G0(j, j) + G0(k + m_hap, k + m_hap) -  G0(i, j) -  G0(i, k + m_hap) - G0(j, k + m_hap);
      D_ijk_higher_order(i,j,k) + lamVals(3, i1) == lamSlack(3, i1) + G0(i, i) + G0(j + m_hap, j + m_hap) + G0(k, k) -  G0(i, j + m_hap) -  G0(i, k) + G0(j + m_hap, k);
      D_ijk_higher_order(i,j,k) + lamVals(4, i1) == lamSlack(4, i1) + G0(i, i) + G0(j + m_hap, j + m_hap) + G0(k + m_hap, k + m_hap) -  G0(i, j + m_hap) -  G0(i, k + m_hap) - G0(j + m_hap, k + m_hap);
      D_ijk_higher_order(i,j,k) + lamVals(5, i1) == lamSlack(5, i1) + G0(i + m_hap, i + m_hap) + G0(j, j) + G0(k, k) -  G0(i + m_hap, j) -  G0(i + m_hap, k) - G0(j, k);
      D_ijk_higher_order(i,j,k) + lamVals(6, i1) == lamSlack(6, i1) + G0(i + m_hap, i + m_hap) + G0(j, j) + G0(k + m_hap, k + m_hap) -  G0(i + m_hap, j) -  G0(i + m_hap, k + m_hap) - G0(j, k + m_hap);
      D_ijk_higher_order(i,j,k) + lamVals(7, i1) == lamSlack(7, i1) + G0(i + m_hap, i + m_hap) + G0(j + m_hap, j + m_hap) + G0(k, k) -  G0(i + m_hap, j + m_hap) -  G0(i + m_hap, k) - G0(j + m_hap, k);
      D_ijk_higher_order(i,j,k) + lamVals(8, i1) == lamSlack(8, i1) + G0(i + m_hap, i + m_hap) + G0(j + m_hap, j + m_hap) + G0(k + m_hap, k + m_hap) -  G0(i + m_hap, j + m_hap) -  G0(i + m_hap, k + m_hap) - G0(j + m_hap, k + m_hap); 
  end;

  minimize val_min;
  subject to
  G0==semidefinite(m);
  lamVals >= 0;
  lamSlack >= 0;
  sum(sum(G0))==0;

  cvx_end


  % Obtain 3D configuration via eigenvector decomposition
  [V,L]=eig(G0);
  x0 = V(:,(m-dim+1):m)*(L((m-dim+1):m, (m-dim+1):m))^0.5;