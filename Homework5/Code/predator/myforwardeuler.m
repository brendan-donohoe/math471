function result = myforwardeuler(fun,tn,yn,h)
  result = yn + h*fun(tn,yn);
end