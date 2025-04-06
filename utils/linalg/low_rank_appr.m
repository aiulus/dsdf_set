function H = low_rank_appr(H)
    [U, S, V] = svd(H);
    rank_cutoff = floor(size(H, 2) * 0.9); % Retain 90% of singular values
    S(rank_cutoff+1:end, rank_cutoff+1:end) = 0;
    H = U * S * V'; % Reconstruct low-rank approximation
end