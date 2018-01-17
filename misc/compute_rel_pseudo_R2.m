function R2 = compute_rel_pseudo_R2(Y, Yhat1, Yhat2)

if(nargin < 4)
  Model1r = sum(Y.*log(eps+Yhat1) - Yhat1);
  Model2r = sum(Y.*log(eps+Yhat2) - Yhat2);
  Sat_r = sum(Y.*log(eps+Y) - Y);


  R2 = (1-(Sat_r-Model2r)./(Sat_r-Model1r))';
end
