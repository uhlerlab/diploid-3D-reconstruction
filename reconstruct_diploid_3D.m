function [cvx_status] = reconstruct_diploid_3D(data_path, rho)


% Read in data - pairwise constraints
pairwise_data = load(strcat(data_path, 'pairwise_distances.mat'));
D_ij_pairwise = pairwise_data.arr;

% Read in homolog-homolog distances
homolog_data = load(strcat(data_path, 'homologue_distances.mat'));
D_ii_homolog = homolog_data.arr;

% Read in data - higher-order constraints
tensor_data = load(strcat(data_path, 'tensor_distances.mat'));
D_ijk_higher_order = tensor_data.arr;

% Read in data - locations of higher-order constraints
constraints_data = load(strcat(data_path, 'constraints.mat'));
indices_ijk_higher_order = constraints_data.arr;

% Read in data - number of domains per chromosome
num_domains_data= load(strcat(data_path, 'num_domains.mat'));
num_domains = num_domains_data.arr;

% Read in data - distances between neighboring beads
inter_domain_data= load(strcat(data_path, 'inter_domain_dist.mat'));
D_neighboring = inter_domain_data.arr;
inter_domain_ind_data= load(strcat(data_path, 'inter_domain_dist_ind.mat'));
indices_neighboring = inter_domain_ind_data.arr;


% Solve SDP
[G0, x0, cvx_status] = solve_sdp(D_ij_pairwise, D_ii_homolog, D_ijk_higher_order, indices_ijk_higher_order, D_neighboring, indices_neighboring, rho)


% Save results
save(strcat(data_path, 'x0.mat'), 'x0');
save(strcat(data_path, 'G0.mat'), 'G0');


% Plot results
m_hap=length(D_ij_pairwise); % total number of homolog domains

plot_3D_coordinates(x0, data_path, m_hap, num_domains)


