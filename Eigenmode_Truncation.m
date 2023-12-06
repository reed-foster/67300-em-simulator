function [Aq, Bq, V, D] = Eigenmode_Truncation(A, B, q)

    [V, D] = eig(A);

    Aq = D((end-q+1):end, (end-q+1):end);
    Bq = V\B;
    Bq = Bq((end-q+1):end);
    
    % q = [2, 5, 10];
    % tic
    % [V, D] = eig(A);o
    % t_eig = toc;
    % 
    % Aq = {};
    % bq = {};
    % cq = {};
    % xq = {};
    % yq = {};
    % 
    % for i = 1:length(q)
    %     tic
    %     % if i == 2
    %     %     profile -memory on %%%%
    %     % end
    %     q(i);
    %     Aq{i} = D((end-q(i)+1):end, (end-q(i)+1):end);
    %     bq{i} = V\B;
    %     bq{i} = bq{i}((end-q(i)+1):end);
    %     cq{i} = cq{i}((end-q(i)+1):end);
    %     % if i == 2
    %     %     profreport %%%%
    %     % end
    %     t_gen(i) = toc;
    % 
    %     xq{i} = Aq{i}\-bq{i};
    %     yq{i} = cq{i}'*xq{i};
    %     syst(i+1) = ss(Aq{i}, bq{i}, cq{i}', 0);
    % end
    % 
    % t_gen_tot = t_gen + t_eig;
end