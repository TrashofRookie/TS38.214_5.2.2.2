function [B, W1, p_WB, beam_max, W2_WB] = beamselect_and_W1form(numberofBeams, N1, N2, Nr, Nt, H_WB)
      if N1 > 1
          O1 = 4;
      else
          O1 = 1;
      end
      if N2 > 1
          O2 = 4;
      else
          O2 = 1;
      end
      CR = zeros(Nr, Nt, O1 * O2);
      response = zeros(O1 * O2, 1);
      D = zeros(Nt, Nt, O1 * O2);
      INDEX = 1;
      for q1 = 0 : (O1 - 1)
          for q2 = 0 : (O2 - 1)
               D(:, :, INDEX) = kron(exp(1i * 2 * pi * (0 : N1-1)' * ((0 : N1-1) * O1 + q1) / (N1 * O1)),exp(1i * 2 * pi * (0 : N2-1)' * ((0 : N2-1) * O2 + q2) / (N2 * O2)));
               CR(:, :, INDEX) = H_WB(:, (1 : Nt)) * D(:, :, INDEX);
               response(INDEX) = max(abs(CR(:, :, INDEX)));
               INDEX = INDEX + 1;
          end
      end
      [~, select1] = max(response);
      b = D(:, :, select1);
      [~, index_b] = maxk(CR(:, :, select1), numberofBeams, 'ComparisonMethod', 'abs');
      B = b(:, index_b);
      W1 = blkdiag(B,B);
      [~, ~, V] = svd(H_WB);
      V_WB = V(:, 1);
      W2_WB = W1' * V_WB;
      [p_WB, beam_max] = wideband_quantization(W2_WB);


end



