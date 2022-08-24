function [B2] = W_beamselect_refine(N1, N2, O1, O2, i_12, L, Nt)
      s_i = 0;

      n = zeros(1,L);
      n_1 = zeros(1,L);
      n_2 = zeros(1,L);
      for i = 0 : (L-1)
          max_x_x = [L-1-i];
          for x_x = L-i : Nt-1-i
              if x_x >= L-i
                  if i_12-s_i >= nchoosek(x_x,L-i)
                      max_x_x = [max_x_x x_x];
                 end
              else
                  if i_12-s_i >= 0
                      max_x_x = [max_x_x x_x];
                  end
              end
          end
          MAX_x_x = max(max_x_x);
          if MAX_x_x >= L-1
              e_i = nchoosek(MAX_x_x, L-1);
          else
              e_i = 0;
          end
          s_i = s_i + e_i;
          n(i + 1) = Nt - 1 - MAX_x_x;
          n_1(i + 1) = mod(N1, n(i + 1));
          n_2(i + 1) = (n(i + 1) - n_1(i + 1))/N1;
      end

%       m_1 = n_1 * O1 + i_11(1);
%       m_2 = n_2 * O2 + i_11(2);

      part1 = repmat((0 : N1-1) * 2 * 1i * pi / (N1 * O1), N2, 1);%指数里的前一项
      part1 = part1(:);
      part2 = repmat((0 : N2-1)' * 2 * 1i * pi / (N2 * O2), 1, N1);%指数里的后一项
      part2 = part2(:);
      B2 = part1 .* (n_1 * O1) + part2 .* (n_2 * O2);

end
% function [W] = W_beamselect(N1,N2,O1,O2,i_11,i_12,L,Nt)
%       s_i = 0;
% 
%       n = zeros(1,L);
%       n_1 = zeros(1,L);
%       n_2 = zeros(1,L);
% %       m_1 = zeros(1,L);
% %       m_2 = zeros(1,L);
%       for i = 0 : (L-1)
%           max_x_x = [L-1-i];
%           for x_x = L-i : Nt-1-i
%               if x_x >= L-i
%                   if i_12-s_i >= nchoosek(x_x,L-i)
%                       max_x_x = [max_x_x x_x];
%                  end
%               else
%                   if i_12-s_i >= 0
%                       max_x_x = [max_x_x x_x];
%                   end
%               end
%           end
%           MAX_x_x = max(max_x_x);
%           if MAX_x_x >= L-1
%               e_i = nchoosek(MAX_x_x, L-1);
%           else
%               e_i = 0;
%           end
%           s_i = s_i + e_i;
%           n(i + 1) = Nt - 1 - MAX_x_x;
%           n_1(i + 1) = mod(N1, n(i + 1));
%           n_2(i + 1) = (n(i + 1) - n_1(i + 1))/N1;
%       end
% %       B = [];
% %       for i = 1 : L
% %           m_1(i) = n_1(i) * O1 + i_11(1);
% %           m_2(i) = n_2(i) * O2 + i_11(2);
% %           B = [B CalcVlm(m_1(i), N1, O1, m_2(i), N2, O2)];
% %       end
%       m_1 = n_1 * O1 + i_11(1);
%       m_2 = n_2 * O2 + i_11(2);
%       B=(0:N1-1)*2i*pi/(N1*O1).*reshape(m_1,1,1,[])+(0:N2-1)'*2i*pi/(N2*O2).*reshape(m_2,1,1,[]);
%       B=reshape(B,[],L);
%       B=exp(B);
% 
%       W=[B,zeros(size(B)); zeros(size(B)), B];
% end