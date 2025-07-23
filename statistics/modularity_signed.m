
function Q = modularity_signed(W, Ci)
    % Computes modularity Q for a signed weighted undirected matrix W
    % and a given partition Ci (cluster labels)
    
    Wp = max(W, 0);   % positive weights
    Wn = -min(W, 0);  % negative weights

    s = sum(W);       % node strength
    s_pos = sum(Wp);
    s_neg = sum(Wn);
    m_pos = sum(s_pos) / 2;
    m_neg = sum(s_neg) / 2;

    Q = 0;
    n = length(Ci);
    
    for i = 1:n
        for j = 1:n
            if Ci(i) == Ci(j)
                % Contribution of positive weights
                Q = Q + (Wp(i,j) - s_pos(i)*s_pos(j)/(2*m_pos));
                % Minus contribution of negative weights
                Q = Q - (Wn(i,j) - s_neg(i)*s_neg(j)/(2*m_neg));
            end
        end
    end
    
    Q = Q / (2*m_pos + 2*m_neg);  % normalize by total weight
end