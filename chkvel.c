void chkvel()
{
  int vort;
  double tmpr2,tmpvx,tmpvy,avgrel,avgabs,tmpdu11,tmpdu12,tmpdu21,
    avgreldu11,avgabsdu11,avgreldu12,avgabsdu12,avgreldu21,avgabsdu21;
   
  avgabs = 0.0;
  avgrel = 0.0;
	
  avgabsdu11 = 0.0;
  avgreldu11 = 0.0;
	
  avgabsdu12 = 0.0;
  avgreldu12 = 0.0;
	
  avgabsdu21 = 0.0;
  avgreldu21 = 0.0;
	
  for (vort=0; vort<N; ++vort)
    {
      tmpr2 = SQR(mblob[vort].blob0.x)+SQR(mblob[vort].blob0.y);
      tmpvx = -0.5*(mblob[vort].blob0.y/tmpr2)*(1.0-exp(-tmpr2/4.0/0.0625));
      tmpvy =  0.5*(mblob[vort].blob0.x/tmpr2)*(1.0-exp(-tmpr2/4.0/0.0625));
      tmpdu11 = mblob[vort].blob0.x*mblob[vort].blob0.y/SQR(tmpr2)*
	(1.0-exp(-tmpr2/4.0/0.0625))-
	0.25*mblob[vort].blob0.x*
	mblob[vort].blob0.y/0.0625/tmpr2*
	exp(-tmpr2/4.0/0.0625);

      tmpdu12 = (-0.5/tmpr2+SQR(mblob[vort].blob0.y/tmpr2))*
	(1.0-exp(-tmpr2/4.0/0.0625))-
	0.25*SQR(mblob[vort].blob0.y)*exp(-tmpr2/4.0/0.0625)/tmpr2/0.0625;

      tmpdu21 = (0.5/tmpr2-SQR(mblob[vort].blob0.x/tmpr2))*
	(1.0-exp(-tmpr2/4.0/0.0625))+
	0.25*SQR(mblob[vort].blob0.x)*exp(-tmpr2/4.0/0.0625)/tmpr2/0.0625;

      if (isnan(tmpvx)) tmpvx = 0.0;
      if (isnan(tmpvy)) tmpvy = 0.0;
	
      if (SQR(mblob[vort].blob0.x)+SQR(mblob[vort].blob0.y)<SQR(0.25))
	{
	  printf("%d %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e  %12.4e\n",
		 vort,
		 sqrt(SQR(tmpvx-mblob[vort].blob0.dx)+
		      SQR(tmpvy-mblob[vort].blob0.dy)),
		 sqrt(SQR(tmpvx-mblob[vort].blob0.dx)+
		      SQR(tmpvy-mblob[vort].blob0.dy))/
		 sqrt(SQR(mblob[vort].blob0.dx)+
		      SQR(mblob[vort].blob0.dy)),
		 sqrt(SQR(tmpdu11-tmpparms[vort].du11)),
		 sqrt(SQR(tmpdu11-tmpparms[vort].du11))/
		 sqrt(SQR(tmpparms[vort].du11)),
		 sqrt(SQR(tmpdu12-tmpparms[vort].du12)),
		 sqrt(SQR(tmpdu12-tmpparms[vort].du12))/
		 sqrt(SQR(tmpparms[vort].du12)),
		 sqrt(SQR(tmpdu21-tmpparms[vort].du21)),
		 sqrt(SQR(tmpdu21-tmpparms[vort].du21))/
		 sqrt(SQR(tmpparms[vort].du21)));
	  avgabsdu11 += sqrt(SQR(tmpdu11-tmpparms[vort].du11));
	  avgreldu11 += sqrt(SQR(tmpdu11-tmpparms[vort].du11))/
	    sqrt(SQR(tmpparms[vort].du11));
	  avgabsdu12 += sqrt(SQR(tmpdu12-tmpparms[vort].du12));
	  avgreldu12 += sqrt(SQR(tmpdu12-tmpparms[vort].du12))/
	    sqrt(SQR(tmpparms[vort].du12));
	  avgabsdu21 += sqrt(SQR(tmpdu21-tmpparms[vort].du21));
	  avgreldu21 += sqrt(SQR(tmpdu21-tmpparms[vort].du21))/
	    sqrt(SQR(tmpparms[vort].du12));
	  avgabs += sqrt(SQR(tmpvx-mblob[vort].blob0.dx)+
			 SQR(tmpvy-mblob[vort].blob0.dy));
	  avgrel += sqrt(SQR(tmpvx-mblob[vort].blob0.dx)+
			 SQR(tmpvy-mblob[vort].blob0.dy))/
	    sqrt(SQR(mblob[vort].blob0.dx)+
		 SQR(mblob[vort].blob0.dy));
	}
    }
  printf("%12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e\n",
	 avgabs/N,avgrel/N,avgabsdu11/N,avgreldu11/N,
	 avgabsdu12/N,avgreldu12/N,
	 avgabsdu21/N,avgreldu21/N);
  exit(1);
}
   
void chkvel2()
{
  int vort;
  double tmpvx,tmpvy,avgrel,avgabs;
   
  avgabs = 0.0;
  avgrel = 0.0;
	
  for (vort=0; vort<N; ++vort)
    {
      dpos_vel(vort);
      tmpvx = mblob[vort].blob0.dx;
      tmpvy = mblob[vort].blob0.dy;
	
      dpos_vel_fast(vort);
	
      if (SQR(mblob[vort].blob0.x)+SQR(mblob[vort].blob0.y)<SQR(0.25))
	{
	  printf("%d %12.4e %12.4e\n",
		 vort,
		 sqrt(SQR(tmpvx-mblob[vort].blob0.dx)+
		      SQR(tmpvy-mblob[vort].blob0.dy)),
		 sqrt(SQR(tmpvx-mblob[vort].blob0.dx)+
		      SQR(tmpvy-mblob[vort].blob0.dy))/
		 sqrt(SQR(mblob[vort].blob0.dx)+
		      SQR(mblob[vort].blob0.dy)));
	  avgabs += sqrt(SQR(tmpvx-mblob[vort].blob0.dx)+
			 SQR(tmpvy-mblob[vort].blob0.dy));
	  avgrel += sqrt(SQR(tmpvx-mblob[vort].blob0.dx)+
			 SQR(tmpvy-mblob[vort].blob0.dy))/
	    sqrt(SQR(mblob[vort].blob0.dx)+
		 SQR(mblob[vort].blob0.dy));
	}
    }
  printf("%12.4e %12.4e\n",avgabs/N,avgrel/N);
  exit(1);
}
   
