function ds2dt = ds2dt(blob,visc)
  ds2dt = visc/2*(blob.a2+1/blob.a2);