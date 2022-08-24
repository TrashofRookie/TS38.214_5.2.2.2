function [B1,B1_polar,p_WB,p_WB_polar,b_quan,b_polar_quan,b_max,b_polar_max] = W1_form_refine(numberofBeams, N1, N2, Nr, Nt, H_WB)
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
      CR_1 = zeros(Nr, Nt, O1 * O2);
      CR_2 = zeros(Nr, Nt, O1 * O2);
      response1 = zeros(O1 * O2, 1);
      response2 = zeros(O1 * O2, 1);
      B = zeros(Nt, Nt, O1 * O2);
      B_polar = zeros(Nt, Nt, O1 * O2);
      INDEX = 1;
      for q1 = 0 : (O1 - 1)
          for q2 = 0 : (O2 - 1)
               B(:,:,INDEX) = kron(exp(1i*2*pi*(0:N1-1)'*((0:N1-1)*O1+q1)/(N1*O1)),exp(1i*2*pi*(0:N2-1)'*((0:N2-1)*O2+q2)/(N2*O2)));
               B_polar(:,:,INDEX) = B(:,:,INDEX);
               CR_1(:, :, INDEX) = H_WB(:, (1 : Nt)) * B(:, :, INDEX);
               CR_2(:, :, INDEX) = H_WB(:, (Nt + 1 : 2 * Nt)) * B_polar(:, :, INDEX);
               response1(INDEX) = max(abs(CR_1(:, :, INDEX)));
               response2(INDEX) = max(abs(CR_2(:, :, INDEX)));
               INDEX = INDEX + 1;
          end
      end
      [~, select1] = max(response1);
      [~, select2] = max(response2);
      b = B(:, :, select1);                                               
      b_polar = B_polar(:, :, select2);
      [p_WB, index_b] = maxk(CR_1(:, :, select1), numberofBeams, 'ComparisonMethod', 'abs');
      [p_WB_polar, index_b_polar] = maxk(CR_2(:, :, select2), numberofBeams, 'ComparisonMethod', 'abs');
      B1 = b(:, index_b);
      B1_polar = b_polar(:, index_b_polar);
      [b_quan,b_max] = wideband_quantization(p_WB);
      [b_polar_quan,b_polar_max] = wideband_quantization(p_WB_polar);
      p_WB = diag(p_WB);
      p_WB_polar = diag(p_WB_polar);


end



