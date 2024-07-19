function est_X = detect_4_PAM(Y,A)

possible_symbols = [-3*A, -A, A, 3*A];

est_X = zeros(size(Y));
for i = 1:length(Y)
    [~, min_distance_symbol] = min(abs(Y(i) - possible_symbols));
    est_X(i) = possible_symbols(min_distance_symbol);
end
end
