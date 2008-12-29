function blob = set_blob(oldblob);

blob=oldblob;

blob.cos2 = cos(blob.th)^2;
blob.sin2 = sin(blob.th)^2;
blob.sincos = sin(blob.th)*cos(blob.th);
blob.sinth = sin(blob.th);
blob.costh = cos(blob.th);
