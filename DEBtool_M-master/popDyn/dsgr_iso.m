function int = dsgr_iso(t, f, rho, rho0, rhoB, lT, lp, li, tp, g, k, vHp, hWG3, hW, hG)
  % modified by Bas Kooijman 2011/04/26, 2015/02/19
  % t: a * kM scaled time
  % int: char equation  S(t) * R(t) * exp(-rt)
  %   assuming that dilution by growth hardly affects surv prob S(t)
  % called by fnsgr_iso

  if hWG3 > 100
    S = exp(- (hW * t).^3 - rho * t);
  else
    hGt = hG * t;
    S = exp(6 * hWG3 * (1 - exp(hGt) + hGt + hGt .^ 2/2) - rho * t);
  end

  l = li - (li - lp) .* exp( - rhoB .* (t - tp));
  int = S .* rho0 .* max(0,(f * l .^ 2 .* (g + lT + l) / (li + g) - k * vHp));
