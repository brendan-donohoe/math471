function result = rungekuttanext(fun,tn,yn,h)
  k1 = fun(tn,yn);
  k2 = fun(tn+(h*0.5),yn+(h*0.5)*k1);
  k3 = fun(tn+(h*0.5),yn+(h.*0.5)*k2);
  k4 = fun(tn+h,yn+h*k3);
  result = yn + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
end