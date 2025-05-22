classdef PAMemoryPoly
    properties
    end

    methods (Access=public)
        function obj = PAMemoryPoly
        end
    end

    methods (Static)
        function coefs = get_coefs_memory_poly(pa_input, pa_output, ployorder, memorydepth)
            x = pa_input(:);
            y = pa_output(:);
            N = numel(x);
            Kstpe = 1;
            K = 1:Kstpe:ployorder;
            Mstep = 1;
            M = 0:Mstep:memorydepth;
            num_K = numel(K);
            num_M = numel(M);
            % get memory polynominal matrix
            x_matrix = PAMemoryPoly.get_maxtirx_memory_ploy(x, K, M);
            max_delay = max(M);
            x_vec = x_matrix(1+max_delay:end,:); % remove delay row
            y_vec = y(1+max_delay:end, :);
            % least square methods
            coefs = x_vec \ y_vec;
        end

        function pa_output = fit_pa_memory_poly(pa_input, coefs, ployorder_vec,  memorydepth_vec)
            x = pa_input;
            N = numel(x);
            K = ployorder_vec;
            M = memorydepth_vec;
            max_delay = max(M);
            x_matrix = PAMemoryPoly.get_maxtirx_memory_ploy(x, K, M);
            if 0
                % remove delay row
                x_matrix(1:max_delay, :) = [];
                y = x_matrix * coefs;
                pa_output = zeros(N,1);
                pa_output(1+max_delay:end,1) = y;
            else
                % non-remove delay row
                y = x_matrix * coefs;
                pa_output = y;
            end
        end

        function x_matrix = get_maxtirx_memory_ploy(x, ployorder_vec, memorydepth_vec)
            N = numel(x);
            K = ployorder_vec;
            M = memorydepth_vec;
            x_matrix = zeros(N, numel(K) * numel(M));
            count = 1;
            for k = K
                for m = M
                    x_m = zeros(N, 1);
                    x_m(1+m:end) = x(1:end-m);
                    abs_x_m_k = abs(x_m).^(k-1);
                    x_matrix(:,count) = x_m .* abs_x_m_k;
                    count = count+1;
                end
            end
        end
    end
end