function [H, H_flat] = custom_hankel(data, order)
    % CUSTOM_HANKEL - Constructs a Hankel matrix with shape 
    %                 [order x (T - order + 1) x m]
    %
    % Inputs:
    %   data    - A T x m matrix, where T is #time steps, and m is the input
    %         dimension.
    %   order   - The desired Hankel matrix order.
    %
    % Output:
    % 
    % H         - A bl. Hankel matrix of size [order x (T - order + 1) x m]
    %
    % H_flat    - Column-wise vertically stacked Hankel matrix of size 
    %            [(order * m) x (T - order + 1)]        

    [m, T] = size(data);

    if order > T
        error("Data matrix has %d rows, but at least %d rows" + ...
            " are required for order %d", T, order, order);
    end
    
    num_columns = T - order + 1; 
    
    % Preallocate Hankel tensor
    H = zeros(order, num_columns, m); 
    
    % Populate the Hankel tensor
    for l=1:order
        %u = data(:, l:(l + num_columns - 1));
        %H(l, :, :) = u;
        for k = 1:num_columns
            H(l, k, :) = data(:, l + k -1);
        end
    end
    
    % Flatten into 2D - Hankel matrix
    H_flat = reshape(permute(H, [1, 3, 2]), order * m, num_columns);
    %disp("Hankel Matrix H of order: "); % DEBUG STATEMENT
    %fprintf("[%d, %d]", size(H, 1), size(H, 2)); % DEBUG STATEMENT
    %disp(H), % DEBUG STATEMENT     
end

