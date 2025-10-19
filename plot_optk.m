af_high_prec = 25 + 5; % Leading term of number of flops in high-precision.
f_mixed_prec = @(ratio,k)(27 / ratio + 8+3*k);
lp_to_hp_ratios = [0.0:0.01:1]; % binary64 to binary16
n_ratios = length(lp_to_hp_ratios);

% Find maximum k such that f_mixed_prec(ratio, k) < f_high_prec.
optk_stit_sylv_orth = floor((19-25*lp_to_hp_ratios + 1-lp_to_hp_ratios)/3);
funk_stit_sylv_orth = (19-25*lp_to_hp_ratios + 1-lp_to_hp_ratios)/3;

optk_ref_sylv_orth = floor((19-25*lp_to_hp_ratios)/3);
funk_ref_sylv_orth = (19-25*lp_to_hp_ratios)/3;

optk_stit_sylv_inv = floor((61/3-25*lp_to_hp_ratios + 1-lp_to_hp_ratios)/3);
funk_stit_sylv_inv = (61/3-25*lp_to_hp_ratios + 1-lp_to_hp_ratios)/3;

optk_ref_sylv_inv = floor((61/3-25*lp_to_hp_ratios)/3);
funk_ref_sylv_inv = (61/3-25*lp_to_hp_ratios)/3;

% Lyapunov equation.
optk_stit_lyap_orth = floor((21-27*lp_to_hp_ratios)/6);
funk_stit_lyap_orth = (21-27*lp_to_hp_ratios)/6;

optk_ref_lyap_orth = floor((19-25*lp_to_hp_ratios)/6);
funk_ref_lyap_orth = (19-25*lp_to_hp_ratios)/6;

optk_stit_lyap_inv = floor((67/3-27*lp_to_hp_ratios)/6);
funk_stit_lyap_inv = (67/3-27*lp_to_hp_ratios)/6;

optk_ref_lyap_inv = floor((61/3-25*lp_to_hp_ratios)/6);
funk_ref_lyap_inv = (61/3-25*lp_to_hp_ratios)/6;

% Plot data.
close
subplot(2,2,1)
hold on
plot(lp_to_hp_ratios, optk_stit_sylv_orth, 'd')
plot(lp_to_hp_ratios, funk_stit_sylv_orth, '-')
plot(lp_to_hp_ratios, optk_ref_sylv_orth, 's')
plot(lp_to_hp_ratios, funk_ref_sylv_orth, '--')
% axis([min(lp_to_hp_ratios), max(lp_to_hp_ratios),...
%       min(optk_stit_inv)-1, max(optk_stit_inv)+1]);
title('Sylvester with re-orthogonalization');
hold off

subplot(2,2,2)
hold on
plot(lp_to_hp_ratios, optk_stit_lyap_orth, 'd')
plot(lp_to_hp_ratios, funk_stit_lyap_orth, '-')
plot(lp_to_hp_ratios, optk_ref_lyap_orth, 's')
plot(lp_to_hp_ratios, funk_ref_lyap_orth, '--')
% axis([min(lp_to_hp_ratios), max(lp_to_hp_ratios),...
%       min(optk_stit_inv)-1, max(optk_stit_inv)+1]);
title('Lyapunov with re-orthogonalization');
hold off

subplot(2,2,3)
hold on
plot(lp_to_hp_ratios, optk_stit_sylv_inv, 'd')
plot(lp_to_hp_ratios, funk_stit_sylv_inv, '-')
plot(lp_to_hp_ratios, optk_ref_sylv_inv, 's')
plot(lp_to_hp_ratios, funk_ref_sylv_inv, '--')
% axis([min(lp_to_hp_ratios), max(lp_to_hp_ratios),...
%       min(optk_stit_inv)-1, max(optk_stit_inv)+1]);
title('Sylvester with inversion');
hold off

subplot(2,2,4)
hold on
plot(lp_to_hp_ratios, optk_stit_lyap_inv, 'd')
plot(lp_to_hp_ratios, funk_stit_lyap_inv, '-')
plot(lp_to_hp_ratios, optk_ref_lyap_inv, 's')
plot(lp_to_hp_ratios, funk_ref_lyap_inv, '--')
% axis([min(lp_to_hp_ratios), max(lp_to_hp_ratios),...
%       min(optk_stit_inv)-1, max(optk_stit_inv)+1]);
title('Lyapunov with inversion');
hold off



% Save data.
outfilename = sprintf('%s/%s', datfolder, 'optk.dat');
outfile = fopen(outfilename, 'w');
header = [' rho ',...
          'optk_stit_sylv_orth    funk_stit_sylv_orth         ',...
          'optk_ref_sylv_orth     funk_ref_sylv_orth    ',...
          'optk_stit_lyap_orth    funk_stit_lyap_orth         ',...
          'optk_ref_lyap_orth     funk_ref_lyap_orth    ',...
          'optk_stit_sylv_inv     funk_stit_sylv_inv         ',...
          'optk_ref_sylv_inv      funk_ref_sylv_inv    ',...
          'optk_stit_lyap_inv     funk_stit_lyap_inv         ',...
          'optk_ref_lyap_inv      funk_ref_lyap_inv\n'];
fprintf(outfile, header);
for i = 1:n_ratios
  fprintf(outfile, ['%.2f ',...
                    '%2d  %.3e    %2d  %.3e    %2d  %.3e    %2d  %.3e    ',...
                    '%2d  %.3e    %2d  %.3e    %2d  %.3e    %2d  %.3e\n'],...
          lp_to_hp_ratios(i),...
          optk_stit_sylv_orth(i), funk_stit_sylv_orth(i),...
          optk_ref_sylv_orth(i), funk_ref_sylv_orth(i),...
          optk_stit_lyap_orth(i), funk_stit_lyap_orth(i),...
          optk_ref_lyap_orth(i), funk_ref_lyap_orth(i),...
          optk_stit_sylv_inv(i), funk_stit_sylv_inv(i),...
          optk_ref_sylv_inv(i), funk_ref_sylv_inv(i),...
          optk_stit_lyap_inv(i), funk_stit_lyap_inv(i),...
          optk_ref_lyap_inv(i), funk_ref_lyap_inv(i));
end
fclose(outfile);
