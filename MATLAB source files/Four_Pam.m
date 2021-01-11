
%% This function calculate the superpositon of m and v.
%% m and v should be two l x 1 vectors. l is the length of the vector.
%% the super position is -3alpha, -alpha, alpha, 3alpha.

function super_position = FourPam_cal_alpha(codewords,alpha)
    m_j1 = codewords(:,1);
    m_j2 = codewords(:,2);
    super_position = [];
    for i = 1:length(m_j1)
        super_position = cat(1,super_position,alpha*2*m_j1(i)+alpha*m_j2(i));        
    end
end