function [h] = plot_3D_coordinates(x0, data_path, m_hap, num_domains)

m=2*m_hap; % total number of domains
dim=3; % dimension
m_chr = length(num_domains); % number of chromsomes

h=figure;
view(3);
set(gca,'GridLineStyle','--');
hold on;

colors = get(gca,'colororder');

if (length(colors) < m_chr)
    colors = distinguishable_colors(m_chr);
end

current_index = 1;
current_index2 = m_hap + 1;

p_list = zeros(m_chr, 1);
legend_chrom = cell(m_chr, 1);

for chrom=1:m_chr
    m_doms = num_domains(chrom);
    
    % get homolog 1
    x1=x0(current_index:(current_index + m_doms-1),1);
    y1=x0(current_index:(current_index + m_doms-1),2);
    z1=x0(current_index:(current_index + m_doms-1),3);
    current_index = current_index + m_doms;

    % get homolog 2
    x2=x0(current_index2:(current_index2 + m_doms-1),1);
    y2=x0(current_index2:(current_index2 + m_doms-1),2);
    z2=x0(current_index2:(current_index2 + m_doms-1),3);

    current_index2 = current_index2 + m_doms;

    % Plot in 3D
    p = plot3(x1,y1,z1,'LineWidth',1.5, 'color', colors(chrom, :));
    p_list(chrom) = p;

    hold on;
    p1 = plot3(x2,y2,z2,'--', 'LineWidth',1.5, 'color', colors(chrom, :));

    legend_chrom{chrom} = int2str(chrom);
end

legend(p_list, legend_chrom);

% Save figure to file
filename = strcat(data_path, 'recostruction_lines.pdf');
print(h, filename,'-dpdf');
