
        program rf4fgood

c       Copyright (C) 2002  Leo Breiman and Adele Cutler

c	This is free open source software but its use in any 
c	commercial product that is sold for profit is prohibited 
c	without the written consent of Leo Breiman and Adele Cutler. 
c	We appreciaate bug notices and suggested improvements.

c	leo@stat.berkeley.edu   adele@math.usu.edu


c       SEE THE GUIDE TO USING RANDOM FORESTS FOR INSTRUCTIONS 


	external scale
c       set all parameters   

c       line #1 information about the data
c       line #2 setting up the run
c       line #3 variable importance options
c       line #4 options using proximities
c	line #5 filling in missing data
c	line #6 setting up parallel coordinates
c	line #7 saving and rerunning the forest
c	line #8 some output controls


	 parameter(
     1	mdim= 27, nsample0= 4, nclass=2,
     &
     1	maxcat=3,ntest= 2,label=1,
     &
     2	jbt= 10, mtry= 10,look=10,ndsize=10,iaddcl=0, jclasswt=1,
     &
     3	imp=0,impfast=1,
     &
     4	ndprox=0,noutlier=0,iscale=0,mdimsc=1,
     &
     5	missquick=0,missright=0,code=-100,
     &
     6	llcoor=0,ncoor=1,
     &
     7	isaverf=1,isavepar=0,irunrf=0,ireadpar=0,
     &
     8	isumout=1,infoutr=0,infouts=1,iproxout=0,iclassout=1,
     &
     &
     &
     &  nsample=(iaddcl+1)*nsample0,
     &	nrnodes=2*(nsample/ndsize)+1,
     &	mimp=imp*(mdim-1)+1,
     &	nimp=imp*(nsample-1)+1,
     &	near=ndprox*(nsample0-1)+1,
     &	nmiss=missright*(nsample-1)+1)
	
	
     
c       MAIN PROGRAM

	
c       DIMENSIONING OF ARRAYS

	
        real x(mdim,nsample),dgini(nrnodes),xbestsplit(nrnodes),
     &	v(nsample),tx(nsample),wl(nclass),classpop(nclass,nrnodes),
     &	rimpmarg(mdim,nsample),tclasscat(nclass,maxcat),fill(mdim),
     &	rmargin(nsample),q(nclass,nsample),win(nsample),tp(nsample),
     &	outlier(near),wr(nclass),diffmarg(mdim),cntmarg(mdim),
     &	rmissimp(mimp),tout(near),countimp(nclass,nimp,mimp),
     &	tmissts(nclass),errimp(mimp),tmiss(nclass),yg(mdim,ncoor),
     &	xc(maxcat),dn(maxcat),cp(maxcat),cm(maxcat),
     &	votecat(maxcat),tclasspop(nclass),xts(mdim,ntest),
     7	qts(nclass,ntest),tgini(mdim),tsgini(mdim),zp(3,mdim),
     &	tw(nrnodes),tn(nrnodes),tnodewt(nrnodes),classwt(nclass),
     &	counttr(nclass,nsample),countts(nclass,ntest)
     	




  	integer cat(mdim),cl(nsample),out(nsample),jests(ntest),
     &	bestsplitnext(nrnodes),treemap(2,nrnodes),bestvar(nrnodes),
     &	bestsplit(nrnodes),nodestatus(nrnodes),nodepop(nrnodes),
     &	parent(nrnodes),jin(nsample),jvr(nsample),mbax(mdim),
     &  nodeclass(nrnodes),nodex(nsample),ncts(nclass),jints(ntest),
     &	nodestart(nrnodes),ta(nsample),ncase(nsample),iv(mdim),
     &	jest(nsample),isort(nsample),ncp(near),clp(near),jts(ntest),
     &	jtr(nsample),nc(nclass),mtab(nclass,nclass),idmove(nsample),
     &	nodexts(ntest),nrcat(maxcat),kcat(maxcat),ncatsplit(maxcat),
     &	clts(ntest)
    
        integer a(mdim,nsample),at(mdim,nsample),b(mdim,nsample)
	integer nbestcat(maxcat,nrnodes),missing(mdim,nmiss)
     	
c     	used in scaling 
      	double precision prox(near,near),y(near),u(near),
     &	dl(mdimsc),xsc(near,mdimsc),red(near)
     
	character*500 text
	character*100 inputfile
	character*100 savedtreefile
	character*100 impfastfile
	character*100 testoutfile
	character*100 testinfile
	
c	read old tree and parameters of original tree run

 	if (irunrf.eq.1)
     &		open(1,file='saveforest',status='old')
 	if (ireadpar.eq.1)
     &		 open(2,file='saveparams', status='old')
	

c	SPECIFY THE NAMES OF OUTPUT FILES
	
c 	ADD by QYJ 
	write(*,*) 'Command: inputfile savetreefile impfastfile testinfile testoutfile'; 
	write(*,*) ''; 
        call getarg(2,savedtreefile)
 	call getarg(3,impfastfile)
	call getarg(5,testoutfile)
	
 	if (isaverf.eq.1)
     &		 open(1, file= savedtreefile )
 	if (isavepar.eq.1)
     &		 open(2, file='saveparams-temp')
 	if (infoutr.eq.1)
     &		 open(3, file='info-temp')
 	if (infouts.eq.1)
     &		open(4, file= testoutfile )
     
 	if (iproxout.eq.1)
     &		open(7, file='prox-temp')
	if (iscale.eq.1)
     &		open(8, file='scale-temp')
 	if (noutlier.eq.1)
     &		open(9, file='outlier-temp')
	if (imp.eq.1)
     &		open(11,file='imp-temp')
 	if (missright.eq.1)
     &		open(13,file='xfilled-temp')
 	if (impfast.eq.1)
     &		open(14,file=impfastfile)
 	if (llcoor.eq.1)
     &		 open(15,file='parall-temp')

c	READ IN DATA

c	read in training data or data to run down a saved forest
	
c	see the manual for format
	
c 	ADD by QYJ 
        call getarg(1,inputfile)

c	read in the data
        open(10,file=inputfile, status='old')
        do k=1,nsample
	       	read(10,*) (x(j,k),j=1,mdim ),cl(k)
        end do
        close(10)

c 	ADD by QYJ
c	read in the test data
	if(ntest.gt.1) then
	
	call getarg(4,testinfile)
        
        open(20,file=testinfile, status='old')
        do k=1,ntest
	       	read(20,*) (xts(j,k),j=1,mdim ),clts(k)
        end do 
        close(20)
 	end if
 
c	see the manual for format

	if(ireadpar.eq.1) goto 888
	if(irunrf.eq.1) goto 999

	write(*,*) 'start'

c       SET CATEGORICAL VALUES+++++++++++++++++++++++++++++++++++
        
        do m=1,mdim
        cat(m)=1
        end do

        if(maxcat.ge.2) then
c       fill in for all variables with cat(m)>1
c        do m= 1, 110
c        	cat(m)=2
c        end do
        end if
	
	if(isaverf.eq.1) then
		write(1,*) '# mdim, nsampleTrain, nclass, nrnodes, jbt, maxcat'
		write(1,*) mdim,nsample,nclass,nrnodes,jbt,maxcat
		write(1,*) (cat(m),m=1,mdim)
	end if
	
c	SET SEED

	do nr=1,344
	zz=rrand()
	end do
	

c	GIVE CLASS WEIGHTS


	do j=1,nclass
	classwt(j)=1
	end do

	
	if(jclasswt.eq.1) then
c	fill in for each class with weight >1
	classwt(1)= 3
	end if
	
	
c	DO PRELIMINARY MISSING DATA
	
	if(missright.eq.1) then
	call zerm(missing,mdim,nsample)
	do n=1,nsample
	do m=1,mdim
	if(x(m,n).eq.code) missing(m,n)=1
	end do
	end do
	nrep=0
	errtag=0
	end if

	if(max(missquick,missright).eq.1) then
	call roughfix(x,v,ncase,mdim,nsample,xts,ntest,cat,code,
     &	nrcat,maxcat,fill)
	if(isaverf.eq.1) write(1,9999) (fill(m), m=1,mdim)
9999	format (1X, F5.2, 1X, F10.6, 1X, F10.6, 8(1X, F5.2), 3(1X, F3.1), F10.6)
	end if
     	
     
c       SET UP DATA TO ADD A CLASS++++++++++++++++++++++++++++++
	
	if(iaddcl.ge.1) then
	call createclass(x,cl,nsample0,nsample,mdim,iaddcl)
	end if
        
	
c	COUNT CLASS POPULATIONS
	
	call zerv(nc,nclass)
	do n=1,nsample
	if(cl(n).lt.1.or.cl(n).gt.nclass) then 
	write(*,*) 'error in class label', n,cl(n)
	stop
	end if
	nc(cl(n))=nc(cl(n))+1
	end do
	if(ntest.gt.1.and.label.eq.1) then
	call zerv(ncts,nclass)
	do n=1,ntest
	ncts(clts(n))=ncts(clts(n))+1
	end do
	end if
	
	if(iclassout.eq.1) then
	write(*,*) 'class counts-training data'
	write(*,*) (nc(j), j=1,nclass)
	if(ntest.gt.1)write(*,*) 'class counts-test data'
	write(*,*) (ncts(j), j=1,nclass)
	end if
	print*
	
	
	
c       INITIALIZE FOR RUN

13	continue
	
 	call zermr(countts,nclass,ntest)
     	call zermr(counttr,nclass,nsample)
	call zerv(out,nsample)
        call zervr(tgini,mdim)
	call zervr(tsgini,mdim)
        call zermd(prox,near,near)
	call zerv(jints,ntest)
	
	
        
        if(imp.eq.1) then
        do m=1,mdim
        do n=1,nsample
        do k=1,nclass
        countimp(k,n,m)=0
        end do
        end do
        end do
        end if
	
	if(jclasswt.eq.0) then
	do k=1,nrnodes
	tnodewt(k)=1
	end do
	end if

	
	call makea(x,mdim,nsample,cat,isort,v,at,b,mbax)
	
C       START RUN   $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	
        do jb=1,jbt
	
	call zerv(nodestatus,nrnodes)
	call zerm(treemap,2,nrnodes)
	call zervr(xbestsplit,nrnodes)
	call zerv(nodeclass,nrnodes)
	call zerv(bestvar,nrnodes)
	
	
        call zerv(jin,nsample)
        call zervr(tclasspop,nclass)
        call zervr(win,nsample)
	
	do n=1,nsample
        k=int(rrand()*nsample)+1
	win(k)=win(k)+classwt(cl(k))
        jin(k)=jin(k)+1
	tclasspop(cl(k))=tclasspop(cl(k))+classwt(cl(k))
	end do
	

	
	
	call eqm(a,at,mdim,nsample)
	
	call moda(a,nuse,nsample,mdim,cat,maxcat,ncase,jin)
	
	
     	call buildtree(a,b,cl,cat,mdim,nsample,nclass,treemap,
     &  bestvar,bestsplit,bestsplitnext,dgini,nodestatus,nodepop,
     &	nodestart,classpop,tclasspop,tclasscat,ta,nrnodes,
     &	idmove,ndsize,ncase,parent,mtry,nodeclass,ndbigtree,
     &	win,wr,wl,nuse,kcat,ncatsplit,xc,dn,cp,cm,maxcat,
     &	nbestcat)
       
	
	call xtranslate(x,mdim,nrnodes,nsample,bestvar,bestsplit,
     &  bestsplitnext,xbestsplit,nodestatus,cat,ndbigtree)
     
     	if(jclasswt.eq.1) then
     	call getweights(x,nsample,mdim,treemap,nodestatus,
     &  xbestsplit,bestvar,nrnodes,ndbigtree,
     &	cat,maxcat,nbestcat,jin,win,tw,tn,tnodewt)
     	end if
     	
     
     	if(isaverf.eq.1) then
	write(1,*) jb,ndbigtree
   	do n=1,ndbigtree
		write(1,*) n,nodestatus(n),bestvar(n), treemap(1,n),
     &		treemap(2,n),nodeclass(n),xbestsplit(n),tnodewt(n), 
     &		(nbestcat(k,n),k=1,maxcat)
	end do
     	end if
	

	
C        GET OUT-OF-BAG ESTIMATES
        
        call testreebag(x,nsample,mdim,treemap,nodestatus,
     &  xbestsplit,bestvar,nodeclass,nrnodes,ndbigtree,
     &	cat,jtr,nodex,maxcat,nbestcat,jin,dgini,tgini)
   
	
	
	
	do n=1,nsample
        if(jin(n).eq.0) then
	counttr(jtr(n),n)=counttr(jtr(n),n)+tnodewt(nodex(n))
        out(n)=out(n)+1
	end if
        end do
	
	
	
c        GET TEST SET ERROR ESTIMATES
        
        
        if(ntest.gt.1) then
        call testreebag(xts,ntest,mdim,treemap,nodestatus,
     &  xbestsplit,bestvar,nodeclass,nrnodes,ndbigtree,
     &	cat,jts,nodexts,maxcat,nbestcat,jints,dgini,tsgini)
    
     	
	do n=1,ntest
	countts(jts(n),n)=countts(jts(n),n)+tnodewt(nodexts(n))
        end do
	end if
	
	
c        DO PROXIMITIES

        if(ndprox.eq.1) then
	do n=1,near
        do k=1,near
        if(nodex(k).eq.nodex(n)) then
	prox(k,n)=prox(k,n)+dble(1)
	end if
        end do
        end do 
        end if
        
     
        
     	
c       DO VARIABLE IMPORTANCE

	if(imp.eq.1) then
        call zerv(iv,mdim)
        do kt=1,ndbigtree
        if(nodestatus(kt).ne.-1) iv(bestvar(kt))=1
        end do
	
        do mr=1,mdim
        if(iv(mr).eq.1) then
	call permobmr(mr,x,tp,tx,jin,nsample,mdim)
	
	call testreebag(x,nsample,mdim,treemap,nodestatus,
     &  xbestsplit,bestvar,nodeclass,nrnodes,ndbigtree,
     &	cat,jvr,nodex,maxcat,nbestcat,jin,dgini,tgini)
    
        do n=1,nsample
	if(jin(n).eq.0)  countimp(jvr(n),n,mr)=
     &	countimp(jvr(n),n,mr)+tnodewt(nodex(n))
        end do
	do n=1,nsample
        x(mr,n)=tx(n)
        end do
	end if
        end do !mr
	do mr=1,mdim
	if(iv(mr).eq.0) then
	do n=1,nsample
	if(jin(n).eq.0) countimp(jtr(n),n,mr)=
     &	countimp(jtr(n),n,mr)+tnodewt(nodex(n))
	end do
	end if
	end do
	end if
	


        
c        GIVE RUNNING OUTPUT

        if(mod(jb,look).eq.0.or.jb.eq.jbt) then
	
	
	call comperrtr(counttr,cl,nsample,nclass,errtr,
     &  tmiss,nc,jest,out)

     	if(iclassout.eq.1) then
     	write(*,'(i8,100f10.2)') 
     &	jb,100*errtr,(100*tmiss(j), j=1,nclass)
     	else
	write(*,'(i8,2f10.2)') jb,100*errtr
	end if
	
	if (ntest.gt.1.and.label.eq.1)then
	call comperrts(countts,clts,ntest,nclass,errts,
     &  tmissts,ncts,jests,label)
	if(iclassout.eq.1) then
	write(*,'(i8,20f10.2)') jb,100*errts,
     &	(100*tmissts(j),j=1,nclass)
     	else
	write(*,'(i8,2f10.2)') jb,100*errts
	end if
     	if(jb.eq.jbt.and.missright.eq.1) errtr0=100*errtr
	end if
	print *
	end if
	
	end do !jb
101	continue
        
c	DEFINE VOTES
	
	do n=1,nsample
	do j=1,nclass
	q(j,n)=counttr(j,n)/out(n)
	end do
	end do
	
	If (ntest.gt.1) then
	do n=1,ntest
	do j=1,nclass
	qts(j,n)=countts(j,n)/jbt
	end do
	end do
	end if

c	FILL IN MISSING VALUES AND TRY AGAIN

	if(missright.eq.1) then
	nrep=nrep+1
	if(nrep.eq.6) stop
	if(errtr0.ge.errtag-.25.and.nrep.ge.2) goto 123
	errtag=errtr0
	write(*,*) 'nrep', nrep
	call impute(x,prox,near,mdim,nsample,nsample0,
     &	missing,maxcat,votecat,cat)
	goto 13
123	continue
	endif

	
	
c       (((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((
	
c       OUTPUT  OUTPUT OUTPUT============================

C       SUMMARY OUTPUT

        if (isumout.eq.1) then
	if(ntest.eq.1.or.label.eq.0)then
        write(*,*) 'final error rate %    ',100*errtr
        else
        write(*,*) 'final error rate %    ',100*errtr
        write(*,*) 'final error test %    ',100*errts
	end if
	print *
        call zerm(mtab,nclass,nclass)
        do n=1,nsample
        if(jest(n).gt.0) mtab(cl(n),jest(n))=mtab(cl(n),jest(n))+1
        end do
        write(*,*) '          true class '
        print *
        write(*,'(20i10)')  (i,i=1,nclass)
        print *
        do j=1,nclass
        write(*,'(20i6)')  j,(mtab(i,j), i=1,nclass)
        end do
        print *
        end if
	
c	PARALLEL COORDINATES
	
	if(llcoor.eq.1) then
	call parcoor(q,v,b,ncase,yg,nsample,mdim,nclass,ncoor,
     &	mbax,zp,cat,x)
     	end if
	
C       SEND INFO ON TRAINING AND/OR TEST SET DATA TO FILE++++++++++++++++++++++++++++ 

        if(infoutr.eq.1) then
        do n=1,nsample
        write(3,'(3i5,50f10.3)') n,cl(n),jest(n),
     &  (q(j,n), j=1,nclass)
        end do
        close(3)
        end if
	
	if(infouts.eq.1) then
	do n=1,ntest
	if (label.eq.1) then
	write(4,'(3i5,50f10.3)') n,clts(n),jests(n),
     &  (qts(j,n), j=1,nclass)
     	else
	write(4,'(3i5,50f10.3)') n,jests(n),
     &  (qts(j,n), j=1,nclass)
     	end if
     	end do
	close(4)
	end if
        
C       SEND PROXIMITY DATA TO FILE++++++++++++++++++++++++++++++++
        
        if(ndprox.eq.1) then
        do k=1,near
        do n=1,near
        prox(n,k)=prox(n,k)/jbt
        if(k.eq.n) prox(n,k)=1
        end do
        end do
	end if
        
        if(iproxout.eq.1) then 
        do n=1,near
        do k=1,near
        write(7,'(2i5,f10.3)') n,k,prox(n,k)
        end do
        end do
        close(5)
        end if
	
C	SEND SCALING DATA TO FILE+++++++++++++++++++++++++++++++++
        
	if(iscale.eq.1) then
	call scale(prox,xsc,y,u,dl,near,mdimsc,red)
	do n=1,mdimsc
        write(*,'(i5,f10.3)') n,dl(n)
	end do
	do n=1,near
        write(8,'(3i5,10f10.3)') n,cl(n),jest(n),(xsc(n,k),k=1,mdimsc)
	end do
	close(8)
	end if

c	SEND IMPUTED DATA TO FILE+++++++++++++++++++++++++++++++++
	
	if(missright.eq.1) then
	if(infoutr.eq.1) then
	do n=1,nsample0
	write(13,*) cl(n), (x(m,n), m=1,mdim)
	end do
	end if
	close(13)
	end if
	
	
C       SEND OUTLIER DATA TO FILE+++++++++++++++++++++++++++++
        
        
        if(noutlier.eq.1) then
	call locateout(prox,cl,near,nsample,nclass,ncp,
     &	iaddcl,outlier,tout,isort,clp)
     	nlier=0
        do n=1,near
	if(outlier(n).ge.10) nlier=nlier+1
	end do
	if (nlier.eq.0) write(*,*) 'no outliers detected'
	if(nlier.gt.0) then
	do n=1,near
	if(outlier(n).ge.10) then
        write(9,'(2i5,f10.3)') ncp(n),clp(n),outlier(n)
     	end if
        end do
	end if
        close(9)
        end if
	
	
	

C	SENDS FASTIMP TO FILE	
	
	avtg=0
	do m=1,mdim
	avtg=avtg+tgini(m)
	end do
	do m=1,mdim
	tgini(m)=100*tgini(m)/avtg
	end do
	
	if(impfast.eq.1) then
	do m=1,mdim
	write(14,*) m,tgini(m)
	end do
	close(14)
	end if
	
c	SENDS STANDAARD IMP OUTPUT TO FILE
        
        if (imp.eq.1) then
	do n=1,nsample
	smax=0
        do j=1,nclass
	if(j.ne.cl(n)) smax=amax1(q(j,n),smax)
  	end do
	pth=q(cl(n),n)
        rmargin(n)=pth-smax
	end do
	
	call finishimp(rmissimp,countimp,out,cl,nclass,mdim,
     &	nsample,errimp,rimpmarg,diffmarg,rmargin)
     	
     	do m=1,mdim
        write(11,'(i5,10f10.2)')m,100*amax1(diffmarg(m),0.0)
        end do
        close(11)
        end if
	


cc	SEND RUN PARAMETERS TO FILE AND READ THEM LATER ++++++++++++++++++++++++++

	if(isavepar.eq.1) then
	write(2,*) nsample0,mdim,maxcat,nclass,jbt,
     &	jclasswt,missquick,missright,code,nrnodes,
     &	100*errtr
     
c		type in comments up to 500 characters long
c		between the ' ' in the line below.

     	 write(2,*) 'this is a test run to verify that my 
     &	 descriptive output works.'
	 close (2)
	 end if

888     continue

	if(ireadpar.eq.1) then

        read(2,*) n0,n1,n2,n3,n4,n5,n6,n7,n8,n9,er
	write(*,*) 'parameters'
	write(*,*) 'nsample=' ,n0
	write(*,*) 'mdim =  ' ,n1
	write(*,*) 'maxcat = ',n2
	write(*,*) 'nclass = ',n3
	write(*,*) 'jbt =   ' ,n4
	write(*,*) 'jclasswt=',n5
	write(*,*) 'missquick=',n6
	write(*,*) 'missright=',n7
        write(*,*) 'code='    ,n8
        write(*,*) 'nrnodes=' ,n9
	print *
	write(*,*) 'out-of-bag error =', er,'%'
	print *

	read(2,'(500a)') text
	write(*,*) text
	stop
        end if

	 
c	RERUN OLD RANDOM FOREST++++++++++++++++++++++++
	
999	continue
	
	if(irunrf.eq.1) then
	
        call zerv(jin,nsample)
	call runforest(mdim,nsample,nclass,maxcat,nrnodes,
     &	label,jbt,cl,x,xbestsplit,counttr,treemap,nbestcat,
     &	nodestatus,cat,nodeclass,jtr,jest,bestvar,nodex,q,
     &	tmiss,nc,fill,missquick,missright,code,
     &	errtr,jin,dgini,tgini,tnodewt) 
	
	if(label.eq.1) then
     	 write(*,'(20f10.2)') 100*errtr,(100*tmiss(j), j=1,nclass)
     	if(infoutr.eq.1) then
        do n=1,nsample
        write(3,'(3i5,50f10.3)') n,cl(n),jest(n),(q(j,n), j=1,nclass)
        end do
        close(3)
        end if
	end if
	end if
	
	end
	

c       END MAIN


c       SUBROUTINE MAKEA

        subroutine makea(x,mdim,nsample,cat,isort,v,a,b,mbax)
        
        real x(mdim,nsample),v(nsample)
        integer cat(mdim),isort(nsample),a(mdim,nsample),
     &	b(mdim,nsample),mbax(mdim)
        
c       submakea constructs the mdim x nsample integer array a.if c	
c	there are less than 32,000 cases, this can be declared integer,
c	otherwise integer*4. For each numerical variable with values x(m,n),
c	n=1,...,nsample, the x-	values are sorted from lowest to highest. 
c	Denote these by xs(m,n).  Then a(m,n) is the case  number in which 
c	xs(m,n) occurs. The b matrix is also contructed here. if the mth 
c	variable is categorical, then a(m,n) is the category of the nth 
c	case number.  

        do 10 mvar=1,mdim

        if (cat(mvar).eq.1) then
                do 20 n=1, nsample
                v(n)=x(mvar,n)
                isort(n)=n
20              continue
                call quicksort(v,isort,1,nsample,nsample)
                
        
        
c       this sorts the v(n) in ascending order. isort(n) is the 
c	case number of that v(n) nth from the lowest (assume the original
c	case numbers are 1,2,...).  
        
		
                do 35 n=1,nsample-1
                n1=isort(n)
                n2=isort(n+1)
                a(mvar,n)=n1
		if(n.eq.1) b(mvar,n1)=1
                if (v(n).lt.v(n+1)) then
                        b(mvar,n2)=b(mvar,n1)+1
                else
                        b(mvar,n2)=b(mvar,n1)
                endif
35              continue
                a(mvar,nsample)=isort(nsample)
		mbax(mvar)=b(mvar,isort(nsample))

        else

                do 40 ncat=1,nsample
                a(mvar,ncat)=nint(x(mvar,ncat))
40              continue

        endif
10      continue

        do 50 n=1,nsample
        isort(n)=n
50      continue

        end
        
        

                   
        
c	SUBROUTINE MODA
	
	subroutine moda(a,nuse,nsample,mdim,cat,maxcat,
     &	ncase,jin)
	integer a(mdim,nsample),cat(mdim),jin(nsample),
     &	ncase(nsample)
	
	nuse=0
	do n=1,nsample
	if(jin(n).ge.1) nuse=nuse+1
	end do
	
	do m=1,mdim
	k=1
	nt=1
	if(cat(m).eq.1) then
	do n=1,nsample
	if(jin(a(m,k)).ge.1) then
	a(m,nt)=a(m,k)
	k=k+1
	else
	do j=1,nsample-k
	if(jin(a(m,k+j)).ge.1) then
	a(m,nt)=a(m,k+j)
	k=k+j+1
	goto 28
	end if
	end do
	end if
28	continue
	nt=nt+1
	if(nt.gt.nuse) goto 37
	end do
37	continue
	end if
	end do
	
	if(maxcat.gt.1) then
	k=1
	nt=1
	do n=1,nsample
	if(jin(k).ge.1) then
	ncase(nt)=k
	k=k+1
	else
	do j=1,nsample-k
	if(jin(k+j).ge.1) then
	ncase(nt)=k+j
	k=k+j+1
	goto 58
	end if
	end do
	end if
58	continue
	nt=nt+1
	if(nt.gt.nuse) goto 85
	end do
85	continue
	end if
	end
        
        

        
c       SUBROUTINE BUILDTREE
        
        subroutine buildtree(a,b,cl,cat,mdim,nsample,nclass,treemap,
     &  bestvar,bestsplit,bestsplitnext,dgini, nodestatus,nodepop,
     &	nodestart,classpop,tclasspop,tclasscat,ta,nrnodes,
     &	idmove,ndsize,ncase,parent,mtry,nodeclass,ndbigtree,
     &	win,wr,wl,nuse,kcat,ncatsplit,xc,dn,cp,cm,maxcat,
     &	nbestcat)
     
     	

c       Buildtree consists of repeated calls to two subroutines, Findbestsplit and Movedata.
c       Findbestsplit does just that--it finds the best split of the current 
c       node.  Movedata moves the data in the split node right and left so that the data
c       corresponding to each child node is contiguous.  

c       The buildtree bookkeeping is different from that in Friedman's original CART program. 
c       ncur is the total number of nodes to date.  nodestatus(k)=1 if the kth node has been split.
c       nodestatus(k)=2 if the node exists but has not yet been split, and =-1 of the node is
c       terminal.  A node is terminal if its size is below a threshold value, or if it is all
c       one class, or if all the x-values are equal.  If the current node k is split, then its
c       children are numbered ncur+1 (left), and ncur+2(right), ncur increases to ncur+2 and
c       the next node to be split is numbered k+1.  When no more nodes can be split, buildtree
c       returns to the main program.

        integer cl(nsample),cat(mdim),ncatsplit(maxcat),
     &  treemap(2,nrnodes),bestvar(nrnodes),nodeclass(nrnodes),
     &	bestsplit(nrnodes), nodestatus(nrnodes),ta(nsample),
     &  nodepop(nrnodes),nodestart(nrnodes),idmove(nsample),
     &	bestsplitnext(nrnodes),ncase(nsample),parent(nrnodes),
     &	kcat(maxcat)
    	
	integer a(mdim,nsample),b(mdim,nsample)
	integer nbestcat(maxcat,nrnodes)
        
        
        real tclasspop(nclass),classpop(nclass,nrnodes),
     &  tclasscat(nclass,maxcat),win(nsample),wr(nclass),
     &	wl(nclass),dgini(nrnodes),xc(maxcat),dn(maxcat),
     &	cp(maxcat),cm(maxcat)
        
        call zerv(nodestatus,nrnodes)
        call zerv(nodestart,nrnodes)
        call zerv(nodepop,nrnodes)
        call zermr(classpop,nclass,nrnodes)
	call zerm(nbestcat,maxcat,nrnodes)
        
      
        
        
        do 20 j=1,nclass
        classpop(j,1)=tclasspop(j)
20      continue

        ncur=1
        nodestart(1)=1
        nodepop(1)=nuse
        nodestatus(1)=2
        
c       start main loop

        do 30 kbuild=1,nrnodes

        if (kbuild.gt.ncur) goto 50
        if (nodestatus(kbuild).ne.2) goto 30
	
	
        
c               initialize for next call to findbestsplit

        ndstart=nodestart(kbuild)
        ndend=ndstart+nodepop(kbuild)-1
        do 40 j=1,nclass
        tclasspop(j)=classpop(j,kbuild)
	
40      continue
        jstat=0
	
	
       
     	call findbestsplit(a,b,cl,mdim,nsample,nclass,cat,
     &  ndstart,ndend,tclasspop,tclasscat,msplit,decsplit,nbest,
     &  ncase,jstat,mtry,win,wr,wl,kcat,ncatsplit,xc,dn,
     &	cp,cm,maxcat)
        
	
        if(jstat.eq.1) then
                nodestatus(kbuild)=-1
                goto 30
        else
                bestvar(kbuild)=msplit
                dgini(kbuild)=decsplit
                if (cat(msplit).eq.1) then
                        bestsplit(kbuild)=a(msplit,nbest)
                        bestsplitnext(kbuild)=a(msplit,nbest+1)
                else
			lcat=cat(msplit)
			do i=1,lcat
			nbestcat(i,kbuild)=ncatsplit(i)
			end do
                endif
        endif
	
c	
c	
        
        
        call movedata(a,ta,mdim,nsample,ndstart,ndend,idmove,ncase,
     &  msplit,cat,nbest,ndendl,ncatsplit,maxcat)
     
        
        
c       leftnode no.= ncur+1, rightnode no. = ncur+2.

        nodepop(ncur+1)=ndendl-ndstart+1
        nodepop(ncur+2)=ndend-ndendl
	nodestart(ncur+1)=ndstart
        nodestart(ncur+2)=ndendl+1
	

c               find class populations in both nodes
        
        do 60 n=ndstart,ndendl
	nc=ncase(n)
	j=cl(nc)
	classpop(j,ncur+1)=classpop(j,ncur+1)+win(nc)
60      continue
        do 70 n=ndendl+1,ndend
	nc=ncase(n)
        j=cl(nc)
        classpop(j,ncur+2)=classpop(j,ncur+2)+win(nc)
70      continue


c       check on nodestatus

        nodestatus(ncur+1)=2
        nodestatus(ncur+2)=2
        if (nodepop(ncur+1).le.ndsize) nodestatus(ncur+1)=-1
        if (nodepop(ncur+2).le.ndsize) nodestatus(ncur+2)=-1
        popt1=0
        popt2=0
        do j=1,nclass
        popt1=popt1+classpop(j,ncur+1)
        popt2=popt2+classpop(j,ncur+2)
        end do
        
        do j=1,nclass
        if (classpop(j,ncur+1).eq.popt1) nodestatus(ncur+1)=-1
        if (classpop(j,ncur+2).eq.popt2) nodestatus(ncur+2)=-1
        end do

        treemap(1,kbuild)=ncur+1
        treemap(2,kbuild)=ncur+2
        parent(ncur+1)=kbuild
        parent(ncur+2)=kbuild
        nodestatus(kbuild)=1
        ncur=ncur+2
        if (ncur.ge.nrnodes) goto 50
        
30      continue
50      continue

        ndbigtree=nrnodes
        do k=nrnodes,1,-1
        if (nodestatus(k).eq.0) ndbigtree=ndbigtree-1
        if (nodestatus(k).eq.2) nodestatus(k)=-1
        end do

        
        do kn=1,ndbigtree
        if(nodestatus(kn).eq.-1) then
        pp=0
        do j=1,nclass
        if(classpop(j,kn).gt.pp) then
        nodeclass(kn)=j
        pp=classpop(j,kn)
        end if
        end do
        end if
        end do
        
        
        
        end


c       SUBROUTINE FINDBESTSPLIT

c       For the best split, msplit is the variable split on. decsplit is the dec. in impurity. 
c       If msplit is numerical, nsplit is the case number of value of msplit split on,
c       and nsplitnext is the case number of the next larger value of msplit.  If msplit is
c       categorical, then nsplit is the coding into an integer of the categories going left.

        subroutine findbestsplit(a,b,cl,mdim,nsample,nclass,cat,
     &  ndstart,ndend,tclasspop,tclasscat,msplit,decsplit,nbest,
     &  ncase,jstat,mtry,win,wr,wl,kcat,ncatsplit,xc,dn,
     &	cp,cm,maxcat)
     	
        integer cl(nsample),cat(mdim), ncase(nsample),
     &	kcat(maxcat),ncatsplit(maxcat),icat(32)
     
     	integer a(mdim,nsample),b(mdim,nsample)
     
     	real tclasspop(nclass),tclasscat(nclass,maxcat),
     &	win(nsample), wr(nclass),wl(nclass),xc(maxcat),
     &	dn(maxcat),cp(maxcat),cm(maxcat)

	external unpack
     	
	
     
c       compute initial values of numerator and denominator of Gini
        
        pno=0
        pdo=0
        do 10 j=1,nclass
        pno=pno+tclasspop(j)*tclasspop(j)
        pdo=pdo+tclasspop(j)
10      continue
        crit0=pno/pdo
        jstat=0
        
	
        
        
c       start main loop through variables to find best split
	
                critmax=-1.0e20
		
		do 20 mv=1,mtry
200             continue
		mvar=int(mdim*rrand())+1
		if(cat(mvar).eq.1) then
                rrn=pno
                rrd=pdo
                rln=0
                rld=0
                call zervr(wl,nclass)
                do j=1,nclass
                wr(j)=tclasspop(j)
		end do
                critvar=-1e20
		
		
	        do 60 nsp=ndstart,ndend-1
                nc=a(mvar,nsp)
	        u=win(nc)
                k=cl(nc)
                rln=rln+u*(2*wl(k)+u)
                rrn=rrn+u*(-2*wr(k)+u)
                rld=rld+u
                rrd=rrd-u
                wl(k)=wl(k)+u
                wr(k)=wr(k)-u
		
		
  		if (b(mvar,nc).lt.b(mvar,a(mvar,nsp+1))) then
		if(k.ne.cl(a(mvar,nsp+1))) then
		if(amin1(rrd,rld).gt.1.0e-5) then
                crit=(rln/rld)+(rrn/rrd)
		if (crit.gt.critmax) then
                                nbest=nsp
                                critmax=crit
				msplit=mvar
			endif
                end if
                end if
		end if
60		continue

 		else
	
c               compute the decrease in impurity given by categorical splits

        	lcat=cat(mvar)
		call zermr(tclasscat,nclass,maxcat)
        	do 70 nsp=ndstart,ndend
                nc=ncase(nsp)
		l=a(mvar,ncase(nsp))
		!write(*,*) l,nc,cl(nc), win(nc),ncase(nsp),mvar
                tclasscat(cl(nc),l)=tclasscat(cl(nc),l)+win(nc)
70              continue
                nnz=0
                do i=1,lcat
                su=0
                do j=1,nclass
                su=su+tclasscat(j,i)
                end do
		dn(i)=su
		nnz=nnz+1
                end do
                if (nnz.eq.1) then
                critmax=-1.0e25
		goto 20
                end if
		
		
		
	
     	if(lcat.gt.10.and.nclass.eq.2 )then
	call catmaxb(tclasscat,tclasspop,xc,cp,cm,kcat,nclass,
     &	lcat,maxcat,ncatsplit,critmax,pdo,nhit,dn)
     	if(nhit.eq.1)then
	msplit=mvar
	end if
	end if
	
	if(lcat.le.10) then
	call catmax(pdo,tclasscat,tclasspop,nclass,lcat,
     &  ncatsp,critmax,nhit,maxcat)
     	if(nhit.eq.1) then
	msplit=mvar
	call unpack(lcat,ncatsp,icat)
	call zerv(ncatsplit,maxcat)
	do k=1,lcat
	ncatsplit(k)=icat(k)
	end do
	end if
	end if
	
	if(lcat.gt.10.and.nclass.gt.2) then
	call catmaxap(tclasscat,tclasspop,kcat,nclass,lcat,
     &	maxcat,ncatsplit,critmax,kbuild,dn,nhit)
     	if(nhit.eq.1) then
	msplit=mvar
	call zerv(ncatsplit,maxcat)
	do l=1,lcat
	if(kcat(l).eq.-1) ncatsplit(l)=1
	end do
     	end if
	end if
	
	end if !cat
	
		
c       this last subroutine returns those categories going left in the best split. 
c       This is coded into a vector of length maxcat (see under subroutine catmax below for details). 
	
20              continue
25		continue                
                decsplit=critmax-crit0
                if (critmax.lt.-1.0e10) jstat=1
		
        end
	
C	SUBROUTINE CATMAXB
	
	subroutine catmaxb(tclasscat,tclasspop,xc,cp,cm,kcat,
     &	nclass,lcat,maxcat,ncatsplit,critmax,pdo,nhit,dn)
     
     	real tclasscat(nclass,maxcat),xc(maxcat),dn(maxcat),
     &	cp(maxcat),cm(maxcat),tclasspop(nclass)
     
     	integer kcat(maxcat),ncatsplit(maxcat)
	
        nhit=0
	
 	do l=1,lcat
	   if(dn(l).gt.0) then
	      xc(l)=tclasscat(1,l)/dn(l)
	   else
	      xc(l)=0
	   end if
	end do
	
	do nk=1,lcat
	   kcat(nk)=nk
	end do
	call quicksort(xc,kcat,1,lcat,maxcat)
	do j=1,nclass
	   cp(j)=0
	   cm(j)=tclasspop(j)
	end do
		
	
        rrd=pdo
        rld=0
        
		
	do n=1,lcat-1
	   rld=rld+dn(kcat(n))
	   rrd=rrd-dn(kcat(n))
	   do k=1,nclass
	      cp(k)=cp(k)+tclasscat(k,kcat(n))
	      cm(k)=cm(k)-tclasscat(k,kcat(n))
	   end do
	   rln=0
	   rrn=0
	   do k=1,nclass
	      rln=rln+cp(k)**2
	      rrn=rrn+cm(k)**2
	   end do
	
	   if (xc(n).lt.xc(n+1)) then
	       if(amin1(rrd,rld).gt.1) then
                 crit=(rln/rld)+(rrn/rrd)
	         if (crit.gt.critmax) then
                    critmax=crit
		    bestsplit=.5*(xc(n)+xc(n+1))
		    nhit=1
		 endif
               end if
            end if
	 end do !n
	
	if (nhit.eq.1) then
	call zerv(ncatsplit,maxcat)
	do l=1,lcat
	if(dn(l).gt.0) then
	xc(l)=tclasscat(1,l)/dn(l)
	else
	xc(l)=0
	end if
	end do
	do i=1,lcat
	if(xc(i).lt.bestsplit) ncatsplit(i)=1
	end do
	end if
	
	
	end

		
        
C       SUBROUTINE CATMAX

   	
        subroutine catmax(pdo,tclasscat,tclasspop,nclass,lcat,
     &  ncatsp,critmax,nhit,maxcat)
        
        
c       this subroutine finds the best categorical split of a categorical variable
c       with lcat categories, nclass classes and tclasscat(j,l) is the number of cases in
c       class j with category value l. The method used is an exhaustive search over all
c       partitions of the category values.  For the two class problem, there is a
c       faster exact algorithm we will add later.  If lcat.ge.10, the exhaustive search
c       gets slow and there is a faster iterative algorithm we can add later.
        
        parameter(jmax=100)
        real tclasscat(nclass,maxcat),tclasspop(nclass),tmpclass(jmax)
        integer icat(32),n

	external unpack
	
        
	nhit=0
        
        do 10 n=1,(2**(lcat-1))-1
        call unpack(lcat,n,icat)
        do j=1,nclass
	tmpclass(j)=0
	do l=1,lcat
	if(icat(l).eq.1) then
        tmpclass(j)=tmpclass(j)+tclasscat(j,l)
        end if
	end do
        end do
        pln=0
        pld=0
        do 40 j=1,nclass
        pln=pln+tmpclass(j)*tmpclass(j)
        pld=pld+tmpclass(j)
40      continue
        prn=0
        do 50 j=1,nclass
        tmpclass(j)=tclasspop(j)-tmpclass(j)
        prn=prn+tmpclass(j)*tmpclass(j)
50      continue
        tdec=(pln/pld)+(prn/(pdo-pld))
        if (tdec.gt.critmax) then
                critmax=tdec
                ncatsp=n
		nhit=1
	endif
10      continue
        
        end
	
		
c	SUBROUTINE CATMAXAP

	subroutine catmaxap(tclasscat,tclasspop,kcat,
     &	nclass,lcat,maxcat,ncatsplit,critmax,kbuild,dn,nhit)
     
     	 parameter(jmax=100)
     
     	real tclasscat(nclass,maxcat),dn(maxcat),ncr0(jmax),
     &	tclasspop(nclass),ncl(jmax),ncr(jmax),ncl0(jmax),ntr,
     &	ntl,ntl0,ntr0
     
     	integer kcat(maxcat),ncatsplit(maxcat)
	
	nhit=0
	
	do l=1,lcat
	if(rrand().le..5) then
	kcat(l)=-1
	else
	kcat(l)=1
	end if
	end do
	
	do j=1,nclass
	ncl(j)=0
	ncr(j)=0
	ntl=0
	ntr=0
	end do
	
	do l=1,lcat
	do j=1,nclass
	if(kcat(l).eq.-1) then
	ncl(j)=ncl(j)+tclasscat(j,l)
	ntl=ntl+dn(l)
	else
	ncr(j)=ncr(j)+tclasscat(j,l)
	ntr=ntr+dn(l)
	end if
	end do
	end do
	
	crit=0
	do j=1,nclass
	crit=crit+(ncl(j)*ncl(j)/ntl)+(ncr(j)*ncr(j)/ntr)
	end do
	
	
	do k=1,1000
	nchange=0
	do l=1,lcat
	do j=1,nclass
	ncl0(j)=ncl(j)
	ncr0(j)=ncr(j)
	ntl0=ntl
	ntr0=ntr
	ncl(j)=ncl(j)+kcat(l)*tclasscat(j,l)
	ncr(j)=ncr(j)-kcat(l)*tclasscat(j,l)
	ntl=ntl+kcat(l)*dn(l)
	ntr=ntr-kcat(l)*dn(l)
	end do
	critnew=0
	do j=1,nclass
	critnew=critnew+(ncl(j)*ncl(j)/ntl)+(ncr(j)*ncr(j)/ntr)
	end do
	
	if(critnew.gt.critmax) then
	critmax=critnew
	nchange=nchange+1
	kcat(l)=-kcat(l)
	nhit=1
	else
	do j=1,nclass
	ncl(j)=ncl0(j)
	ncr(j)=ncr0(j)
	ntl=ntl0
	ntr=ntr0
	end do
	end if
	end do !l
	if(nchange.eq.0) goto 101
	end do !k
	
101	continue
	
	
	end
	


        
c       SUBROUTINE MOVEDATA     

c       This subroutine is the heart of the buildtree construction. Based on the best split 
c       the data in the part of the a matrix corresponding to the current node is moved to the 
c       left if it belongs to the left child and right if it belongs to the right child.

        subroutine movedata(a,ta,mdim,nsample,ndstart,ndend,idmove,
     &  ncase,msplit,cat,nbest,ndendl,ncatsplit,maxcat)
        integer a(mdim,nsample),ta(nsample),idmove(nsample),
     &  ncase(nsample),cat(mdim),ncatsplit(maxcat)
        
c       compute idmove=indicator of case nos. going left

                
        if (cat(msplit).eq.1) then
                do 10 nsp=ndstart,nbest
                nc=a(msplit,nsp)
                idmove(nc)=1
10              continue
                do 20 nsp=nbest+1, ndend
                nc=a(msplit,nsp)
                idmove(nc)=0
20              continue
                ndendl=nbest
        else
                ndendl=ndstart-1
                do 30 nsp=ndstart,ndend
		nc=ncase(nsp)
                if (ncatsplit(a(msplit,nc)).eq.1) then
                        idmove(nc)=1
                        ndendl=ndendl+1
                else
                        idmove(nc)=0
                endif
30              continue
        endif
        
c       shift case. nos. right and left for numerical variables.        
        
        do 40 msh=1,mdim
        if (cat(msh).eq.1) then
                k=ndstart-1
                do 50 n=ndstart,ndend
                ih=a(msh,n)
                if (idmove(ih).eq.1) then
                        k=k+1
                        ta(k)=a(msh,n)
                endif
50      continue
        do 60 n=ndstart,ndend
        ih=a(msh,n)
        if (idmove(ih).eq.0) then 
                k=k+1
                ta(k)=a(msh,n)
        endif
60      continue

        do 70 k=ndstart,ndend
        a(msh,k)=ta(k)
70      continue
        endif

40      continue
	ndo=0
	if(ndo.eq.1) then
	 do 140 msh=1,mdim
        if (cat(msh).gt.1) then
                k=ndstart-1
                do 150 n=ndstart,ndend
                ih=ncase(n)
                if (idmove(ih).eq.1) then
                        k=k+1
                        ta(k)=a(msh,ih)
                endif
150      continue
        do 160 n=ndstart,ndend
        ih=ncase(n)
        if (idmove(ih).eq.0) then 
                k=k+1
                ta(k)=a(msh,ih)
        endif
160      continue

        do 170 k=ndstart,ndend
        a(msh,k)=ta(k)
170      continue
        endif

140      continue
	end if

	

c       compute case nos. for right and left nodes.
        
        if (cat(msplit).eq.1) then
                do 80 n=ndstart,ndend
                ncase(n)=a(msplit,n)
80              continue
        else
                k=ndstart-1
                do 90 n=ndstart, ndend
                if (idmove(ncase(n)).eq.1) then
                        k=k+1
                        ta(k)=ncase(n)
                endif
90              continue
                do 100 n=ndstart,ndend
                if (idmove(ncase(n)).eq.0) then
                        k=k+1
                        ta(k)=ncase(n)
                endif
100             continue
                do 110 k=ndstart,ndend
                ncase(k)=ta(k)
110             continue
        endif
                
        end
        

c       SUBROUTINE XTRANSLATE

c       this subroutine takes the splits on numerical variables and translates them
c       back into x-values.  It also unpacks each categorical split into a 32-dimensional 
c       vector with components of zero or one--a one indicates that the corresponding
c       category goes left in the split.

        subroutine xtranslate(x,mdim,nrnodes,nsample,bestvar,
     &  bestsplit,bestsplitnext,xbestsplit,nodestatus,cat,
     &	ndbigtree)
    
        integer cat(mdim),bestvar(nrnodes),bestsplitnext(nrnodes),
     &	nodestatus(nrnodes),bestsplit(nrnodes)
        real x(mdim,nsample),xbestsplit(nrnodes)

        do 10 k=1,ndbigtree
        if (nodestatus(k).eq.1) then
                m=bestvar(k)
                if (cat(m).eq.1) then
                        xbestsplit(k)=(x(m,bestsplit(k))+
     &                                 x(m,bestsplitnext(k)))/2
     		else
			xbestsplit(k)=real(bestsplit(k))
                endif
        endif
10      continue
        end
	
c	SUBROUTINE GETWEIGHTS
	
	subroutine getweights(x,nsample,mdim,treemap,nodestatus,
     &  xbestsplit,bestvar,nrnodes,ndbigtree,
     &	cat,maxcat,nbestcat,jin,win,tw,tn,tnodewt)
   
        real x(mdim,nsample),xbestsplit(nrnodes),
     &	win(nsample),tw(nrnodes),tn(nrnodes),tnodewt(nrnodes)
   
	
        integer treemap(2,nrnodes),bestvar(nrnodes),
     &  cat(mdim),nodestatus(nrnodes),jin(nsample)
     
     	integer nbestcat(maxcat,nrnodes)
	
	call zervr(tw,nrnodes)
	call zervr(tn,nrnodes)
	
        do n=1,nsample
	if(jin(n).ge.1) then
        kt=1
        do k=1,ndbigtree
	if (nodestatus(kt).eq.-1) then
	tw(kt)=tw(kt)+win(n)
	tn(kt)=tn(kt)+jin(n)
	goto 100
        end if
	
        m=bestvar(kt)
	
	if (cat(m).eq.1) then
                if (x(m,n).le.xbestsplit(kt)) then 
                        kt=treemap(1,kt)
			else
                        kt=treemap(2,kt)
			endif
		
        else
                jcat=nint(x(m,n))
		
                if (nbestcat(jcat,kt).eq.1) then
                        kt=treemap(1,kt)
                else
                        kt=treemap(2,kt)
                endif
	 endif
	 
	end do
100     continue
	end if
        end do
	do n=1,nrnodes
	if(nodestatus(n).eq.-1) tnodewt(n)=tw(n)/tn(n)
	end do
	end
	
	subroutine testreebag(x,nts,mdim,treemap,nodestatus,
     &  xbestsplit,bestvar,nodeclass,nrnodes,ndbigtree,
     &	cat,jtr,nodex,maxcat,nbestcat,jin,dgini,tgini)
    
        real x(mdim,nts),xbestsplit(nrnodes),dgini(nrnodes),
     &	tgini(mdim)
	
        integer treemap(2,nrnodes),bestvar(nrnodes),
     &  nodeclass(nrnodes),cat(mdim),nodestatus(nrnodes),
     &  nodex(nts),jin(nts),jtr(nts)
     
     	integer nbestcat(maxcat,nrnodes)
        
        call zerv(jtr,nts)
        call zerv(nodex,nts)
	
        do n=1,nts

        kt=1
        do k=1,ndbigtree
	   if (nodestatus(kt).eq.-1) then
	      jtr(n)=nodeclass(kt)
      	      nodex(n)=kt
	      goto 100
           end if
	
           m=bestvar(kt)
	
	   if(jin(n).eq.0) tgini(m)=tgini(m)+dgini(kt)

           if (cat(m).eq.1) then
                if (x(m,n).le.xbestsplit(kt)) then 
                        kt=treemap(1,kt)
		else
                        kt=treemap(2,kt)
		endif
		
           else
                jcat=nint(x(m,n))
		
                if (nbestcat(jcat,kt).eq.1) then
                        kt=treemap(1,kt)
                else
                        kt=treemap(2,kt)
                endif
	   endif
        end do
100     continue
        end do
	end
        
C       SUBROUTINE COMPERRTR
        
        subroutine comperrtr(counttr,cl,nsample,nclass,errtr,
     &  tmiss,nc,jest,out)
    
        integer cl(nsample),nc(nclass),jest(nsample),out(nsample)
	real counttr(nclass,nsample),tmiss(nclass)
      
        call zervr(tmiss,nclass)
	errtr=0
	
	do n=1,nsample
	cmax=0
        do j=1,nclass
	ctemp=counttr(j,n)/(out(n)+.001)
        if (ctemp.gt.cmax) then
        jmax=j
        cmax=ctemp
        end if
        end do
	jest(n)=jmax
	if (jmax.ne.cl(n)) then
	tmiss(cl(n))=tmiss(cl(n))+1
	errtr=errtr+1
	end if
	end do
	errtr=errtr/nsample
	do j=1,nclass
	tmiss(j)=tmiss(j)/nc(j)
	 end do
	
        end
	
c	SUBROUTINE  COMPERRTS 

	subroutine comperrts(countts,clts,ntest,nclass,errts,
     &  tmissts,ncts,jests,label)
    
        integer clts(ntest),ncts(nclass),jests(ntest)
	real countts(nclass,ntest),tmissts(nclass)
      
        call zervr(tmissts,nclass)
	errts=0
	
	do n=1,ntest
	cmax=0
        do j=1,nclass
	if (countts(j,n).gt.cmax) then
        jmax=j
        cmax=countts(j,n)
        end if
        end do
	jests(n)=jmax
	if (label.eq.1) then
	if (jmax.ne.clts(n)) then
	tmissts(clts(n))=tmissts(clts(n))+1
	errts=errts+1
	end if
	end if
	end do
        errts=errts/ntest
	do j=1,nclass
	tmissts(j)=tmissts(j)/ncts(j)
	 end do
	
	 end
	
	
	
	
C	SUBROUTINE CREATECLASS

	subroutine createclass(x,cl,ns,nsample,mdim,iaddcl)
	real x(mdim,nsample)
     	integer cl(nsample)
	
	do n=1,ns
	cl(n)=1
	end do
	do n=ns+1,nsample
	cl(n)=2
	end do
	 
	if(iaddcl.eq.1) then
	do n=ns+1,nsample
	do m=1,mdim
	k=int(rrand()*ns)+1
	x(m,n)=x(m,k)
	end do
	end do
	end if
	
	end
	


C       SUBROUTINE PERMOBAR
        
        subroutine permobmr(mr,x,tp,tx,jin,nsample,mdim)
        real x(mdim,nsample), tp(nsample),tx(nsample)
	integer jin(nsample)
        kout=0
        call zervr(tp,nsample)
        do n=1,nsample
        if(jin(n).eq.0) then
        kout=kout+1
        tp(kout)=x(mr,n)
        end if
        end do 
        call perm1(kout,nsample,tp)
        iout=0
        do n=1,nsample
        tx(n)=x(mr,n)
        if(jin(n).eq.0) then
        iout=iout+1
        x(mr,n)=tp(iout)
        end if
        end do
        end
        
	
C       SUBROUTINE FINISHIMP
        
        subroutine finishimp(rmissimp,countimp,out,cl,nclass,mdim,
     &  nsample, errimp,rimpmarg,diffmarg,rmargin)
     
        integer cl(nsample),out(nsample)
     
        real errimp(mdim),rimpmarg(mdim,nsample),rmissimp(mdim),
     &  diffmarg(mdim),rmargin(nsample),countimp(nclass,nsample,mdim)
     
     
    	call zervr(rmissimp,mdim)
        call zervr(errimp,mdim)
	do n=1,nsample
	do m1=1,mdim
	clsum=0
	do j=1,nclass
	clsum=clsum+(countimp(j,n,m1)/out(n))
	end do
	do j=1,nclass
	countimp(j,n,m1)=countimp(j,n,m1)/(clsum*out(n))
	end do
	end do
	end do
	
        do m1=1,mdim
	do n=1,nsample
        rmax=0
        smax=0
        do j=1,nclass
        if (countimp(j,n,m1).gt.smax) then
        smax=countimp(j,n,m1)
        imax=j
        end if
	if(j.ne.cl(n)) rmax=amax1(rmax,countimp(j,n,m1))
        end do
	if(imax.ne.cl(n)) rmissimp(m1)=rmissimp(m1)+1
	rimpmarg(m1,n)=countimp(cl(n),n,m1)-rmax
	end do 
        end do 
        do m1=1,mdim
        errimp(m1)=real(rmissimp(m1))/nsample
	end do
        do m=1,mdim
        diffmarg(m)=0
        do n=1,nsample
        diffmarg(m)=diffmarg(m)+(rmargin(n)-rimpmarg(m,n))
        end do
	diffmarg(m)=diffmarg(m)/nsample
	end do
	end
        
C	SUBROUTINE PARCOOR	

	subroutine parcoor(q,v,b,ncase,yg,nsample,mdim,
     &	nclass,ncoor,mbax,zp,cat,x)
	real q(nclass,nsample),v(nsample),yg(mdim,ncoor),
     &	zp(3,mdim),x(mdim,nsample)
     	integer ncase(nsample),b(mdim,nsample),mbax(mdim),
     &	cat(mdim)
	
	
	
	do j=1,nclass
	do n=1,nsample
	v(n)=q(j,n)
	end do
	call quicksort(v,ncase,1,nsample,nsample)
	th=v(nsample-ncoor+1)
	nt=0
	do n=1,nsample
	if(q(j,n).ge.th) then
	nt=nt+1
	ncase(nt)=n
	end if
	end do
	
	do m=1,mdim
	if(cat(m).eq.1) then
	do n=1,min(nt,ncoor)
	yg(m,n)=real(b(m,ncase(n)))/mbax(m)
	end do
	else
	ncat=cat(m)
	do n=1,min(nt,ncoor)
	yg(m,n)=(x(m,ncase(n))-1)/(ncat-1)
	end do
	end if
	end do
	
	do m=1,mdim
	if(cat(m).eq.1) then
	do n=1,ncoor
	v(n)=yg(m,n)
	end do
	call quicksort(v,ncase,1,ncoor,nsample)
	zp(1,m)=v(nint(.25*ncoor))
	zp(2,m)=v(nint(.5*ncoor))
	zp(3,m)=v(nint(.75*ncoor))
	else
	do l=1,3
	zp(l,m)=0
	end do
	end if
	end do
	
	
     	do m=1,mdim
	write(15,*) j,m,(yg(m,n),n=1,ncoor),(zp(i,m),i=1,3)
	end do
	end do
	close(15)

	
	end
	
C	SUBROUTINE LOCATEOUT
 	
	subroutine locateout(prox,cl,near,nsample,nclass,ncp,
     &	iaddcl,outlier,tout,isort,clp)
     
	double precision prox(near,near)
	real outlier(near),tout(near)
	integer ncp(near),cl(nsample),isort(nsample),ntt(0:30),
     &	clp(near)
	
	call zervr(outlier,near)
	
	if (iaddcl.eq.1) then
	jpclass=1
	else
	jpclass=nclass
	end if
	ntt(0)=0
	nt=0
	do jp=1,jpclass
	
	do n=1,near
	if(cl(n).eq.jp) then
	nt=nt+1
	ncp(nt)=n
	end if
	end do
	ntt(jp)=nt
	n1=ntt(jp-1)+1
	n2=ntt(jp)
	
	do i=n1,n2
	outlier(i)=0
	do j=n1,n2
	if(j.ne.i) then
	outlier(i)=outlier(i)+(real(prox(ncp(i),ncp(j))))**2
	end if
	end do
	end do
	
	do i=n1,n2
        if (outlier(i).gt.0) then
	outlier(i)=1.0/outlier(i)
        else
        outlier(i)=1000
        end if
	tout(i)=outlier(i)
	clp(i)=jp
	end do
	
	call quicksort(tout,isort,n1,n2,nsample)
	rmed=tout((n1+n2)/2)
	dev=0
	do i=n1,n2
	dev=dev+abs(tout(i)-rmed)
	end do
	dev=dev/(n2-n1+1)
	do i=n1,n2
	outlier(i)=(outlier(i)-rmed)/dev
	outlier(i)=amax1(outlier(i),0.0)
	end do
	
	end do
	
	end
	
	
	 

c       SUNROUTINE SCALE

	subroutine scale(s,xsc,y,u,dl,ns,nn,red)
	implicit double precision (a-h,o-z)
	double precision s(ns,ns),y(ns),u(ns),dl(ns),
     &	xsc(ns,nn),red(ns)	
	
	do j=1,ns
	red(j)=0
	do i=1,ns
	red(j)=red(j)+s(i,j)
	end do
	red(j)=red(j)/ns
	end do
	
	sred=0
	do j=1,ns
	sred=sred+red(j)
	end do
	sred=sred/ns
	
	do j=1,ns
	do i=1,ns
	s(i,j)=(s(i,j)-red(j)-red(i)+sred)/2
	end do
	end do
	
      
      	do it=1,nn
	
	do n=1,ns
	if(mod(n,2).eq.0) then
	y(n)=1
	else
	y(n)=-1
	end if
	end do
	
	do jit=1,1000
	y2=0
	do n=1,ns
	y2=y2+y(n)**2
	end do
	y2=dsqrt(y2)
	do n=1,ns
	u(n)=y(n)/y2
	end do
	
	
	do n=1,ns
	y(n)=0
	do k=1,ns
	y(n)=y(n)+s(n,k)*u(k)
	end do
	end do
	
	ra=0
	do n=1,ns
	ra=ra+y(n)*u(n)
	end do 
	
	ynorm=0
	do n=1,ns
	ynorm=ynorm+(y(n)-ra*u(n))**2
	end do
	if(ynorm.lt.ra*	1.0e-7)then
	do n=1,ns
	xsc(n,it)=dsqrt(ra)*u(n)
	end do
	dl(it)=ra
	goto 101
	end if
	end do
101	continue
	do k=1,ns
	do n=1,ns
	s(n,k)=s(n,k)-ra*u(n)*u(k)
	end do
	end do
	end do !nn
	end
	
	subroutine roughfix(x,v,ncase,mdim,nsample,
     &	xts,ntest,cat,code,nrcat,maxcat,fill)
     
     	real x(mdim,nsample), v(nsample),fill(mdim)
     	real xts(mdim,ntest)
	integer ncase(nsample),cat(mdim),nrcat(maxcat)
	
	
	do m=1,mdim
	    if(cat(m).eq.1) then
	       nt=0
	       do n=1,nsample
	          if(x(m,n).ne.code) then
	             nt=nt+1
	             v(nt)=x(m,n)
	          end if
	       end do
	       call quicksort (v,ncase,1,nt,nsample)
	       if(nt.gt.0) then
	          rmed=v(nt/2)
	          fill(m)=rmed
	       else
	          rmed=0
	       end if
	       do n=1,nsample
	          if(x(m,n).eq.code) x(m,n)=rmed
	       end do
               if(ntest.gt.1) then
	          do n=1,ntest
	             if(xts(m,n).eq.code) xts(m,n)=rmed
	          end do
	       end if
	    end if
	
	    if (cat(m).gt.1) then
	       lcat=cat(m)
	       call zerv(nrcat,maxcat)
	       do n=1,nsample
	          if(x(m,n).ne.code) then
	             j=nint(x(m,n))
	             nrcat(j)=nrcat(j)+1
	          end if
	       end do
	       nmax=0
	       jmax=1
	       do j=1,lcat
	          if(nrcat(j).gt.nmax) then
	             nmax=nrcat(j)
	             jmax=j
	          end if
	       end do
	       fill(m)=real(jmax)
	       do n=1,nsample
	          if(x(m,n).eq.code) x(m,n)=real(jmax)
	       end do
               if(ntest.gt.1) then
	          do n=1,ntest
	             if(xts(m,n).eq.code) xts(m,n)=fill(m)
	          end do
	       end if
            endif
	 
	 end do !m
	 end
	 


c	SUBROUTINE IMPUTE
	 
	 subroutine impute(x,prox,near,mdim,nsample,
     &	 nsample0,missing,maxcat,votecat,cat)
     	
	 real votecat(maxcat),x(mdim,nsample)
     	 double precision  prox(near,near)
	 integer missing(mdim,nsample),cat(mdim)
     
	 
	 
	 do m=1,mdim
	 do n=1,nsample0
	 if(cat(m).eq.1) then
	 if(missing(m,n).eq.1) then
	 sx=0
	 dt=0
	 do k=1,nsample0
	 if (missing(m,k).ne.1) then
	 sx=sx+real(prox(n,k))*x(m,k)
	 dt=dt+real(prox(n,k))
	 end if
	 end do
	 
	 if(dt.gt.0) x(m,n)=sx/dt
	 end if
	 end if
	 
	 
	 if(cat(m).gt.1) then
	 if(missing(m,n).eq.1) then
	 call zervr(votecat,maxcat)
	 do k=1,nsample0
	 if (missing(m,k).ne.1) then
	 j=nint(x(m,k))
	 votecat(j)=votecat(j)+real(prox(n,k))
	 end if
	 end do !k
	 rmax=0
	 do i=1,cat(m)
	 if(votecat(i).gt.rmax) then
	 rmax=votecat(i)
	 jmax=i
	 end if
	 end do
	 x(m,n)=real(jmax)
	 end if
	 end if
	 end do !n
	 end do !m
	
	 end
	 


	
C	SUBROUTINE RUNFOREST

c	reads a forest file and runs new data through it
	
	subroutine runforest(mdim,nsample,nclass,maxcat,nrnodes,
     &	label,jbt,cl,x,xbestsplit,counttr,treemap,nbestcat,
     &	nodestatus,cat,nodeclass,jtr,jest,bestvar,nodex,q,
     &	tmiss,nc,fill,missquick,missright,code,errtr,jin,dgini,
     &	tgini,tnodewt)
    
	
	
	real x(mdim,nsample),xbestsplit(nrnodes),tmiss(nclass),
     &	q(nclass,nsample),fill(mdim),dgini(nrnodes),
     &	tgini(mdim),counttr(nclass,nsample),
     &	tnodewt(nrnodes)
     
     
	integer treemap(2,nrnodes),nodestatus(nrnodes),
     &	cat(mdim),nodeclass(nrnodes),jest(nsample),
     &	bestvar(nrnodes),jtr(nsample),cl(nsample),
     &	nodex(nsample),nc(nclass),jin(nsample)
     
    	integer nbestcat(maxcat,nrnodes)
     	
	
	
	call zermr(counttr,nclass,nsample)
	
	read(1,*) (cat(m),m=1,mdim)
	
	if(max(missquick,missright).eq.1) then
	read(1,*) (fill(m),m=1,mdim)
	do n=1,nsample
	do m=1,mdim
	if(x(m,n).eq.code) x(m,n)=fill(m)
	end do
	end do
	end if
	
	if(label.eq.1) then
	do n=1,nsample
	nc(cl(n))=nc(cl(n))+1
	end do
	end if
	
	do jb=1,jbt
	read(1,*) ndbigtree
	do n=1,ndbigtree
	read(1,*) k,nodestatus(n),bestvar(n),
     &	treemap(1,n),treemap(2,n),nodeclass(n),
     &	xbestsplit(n),tnodewt(n), (nbestcat(k,n),k=1,maxcat)
     	end do
	
	call testreebag(x,nsample,mdim,treemap,nodestatus,
     &  xbestsplit,bestvar,nodeclass,nrnodes,ndbigtree,
     &	cat,jtr,nodex,maxcat,nbestcat,jin,dgini,tgini)
    
	
	do n=1,nsample
	counttr(jtr(n),n)=counttr(jtr(n),n)+tnodewt(nodex(n))
	end do

	call comperrts(counttr,cl,nsample,nclass,errtr,
     &  tmiss,nc,jest,label)

     	end do
	
     	do n=1,nsample
	qsum=0
	do j=1,nclass
	qsum=qsum+counttr(j,n)
	end do
	do j=1,nclass
	q(j,n)=counttr(j,n)/qsum
	end do
	end do
	
	end
	
	

	
	

C       SUBROUTINE QUICKSORT
        
        subroutine quicksort (v,iperm,ii,jj,kk)
c
c     puts into iperm the permutation vector which sorts v into
c     increasing order.  only elementest from ii to jj are considered.
c     array iu(k) and array il(k) permit sorting up to 2**(k+1)-1 elements
c
c     this is a modification of acm algorithm #347 by r. c. singleton,
c     which is a modified hoare quicksort.
c
      real v(kk)
      integer t,tt,iperm(kk),iu(32),il(32)
     
      m=1
      i=ii
      j=jj
 10   if (i.ge.j) go to 80
 20   k=i
      ij=(j+i)/2
      t=iperm(ij)
      vt=v(ij)
      if (v(i).le.vt) go to 30
      iperm(ij)=iperm(i)
      iperm(i)=t
      t=iperm(ij)
      v(ij)=v(i)
      v(i)=vt
      vt=v(ij)
 30   l=j
      if (v(j).ge.vt) go to 50
      iperm(ij)=iperm(j)
      iperm(j)=t
      t=iperm(ij)
      v(ij)=v(j)
      v(j)=vt
      vt=v(ij)
      if (v(i).le.vt) go to 50
      iperm(ij)=iperm(i)
      iperm(i)=t
      t=iperm(ij)
      v(ij)=v(i)
      v(i)=vt
      vt=v(ij)
      go to 50
 40   iperm(l)=iperm(k)
      iperm(k)=tt
      v(l)=v(k)
      v(k)=vtt
 50   l=l-1
      if (v(l).gt.vt) go to 50
      tt=iperm(l)
      vtt=v(l)
 60   k=k+1
      if (v(k).lt.vt) go to 60
      if (k.le.l) go to 40
      if (l-i.le.j-k) go to 70
      il(m)=i
      iu(m)=l
      i=k
      m=m+1
      go to 90
 70   il(m)=k
      iu(m)=j
      j=l
      m=m+1
      go to 90
 80   m=m-1
      if (m.eq.0) return
      i=il(m)
      j=iu(m)
 90   if (j-i.gt.10) go to 20
      if (i.eq.ii) go to 10
      i=i-1
 100  i=i+1
      if (i.eq.j) go to 80
      t=iperm(i+1)
      vt=v(i+1)
      if (v(i).le.vt) go to 100
      k=i
 110  iperm(k+1)=iperm(k)
      v(k+1)=v(k)
      k=k-1
      if (vt.lt.v(k)) go to 110
      iperm(k+1)=t
      v(k+1)=vt
      go to 100
      end
  	

c       MISCELLANOUS SMALL SUBROUTINES

	 subroutine unpack(l,npack,icat)
        integer icat(32),npack
        do j=1,32
        icat(j)=0
        end do
        n=npack
        icat(1)=mod(n,2)
        do 10 k=2,l
        n=(n-icat(k-1))/2
        icat(k)=mod(n,2)
10      continue
        end

	 subroutine zerv(ix,m1)
        integer ix(m1)
        do 10 n=1,m1
        ix(n)=0
10      continue
        end
        
        subroutine zervr(rx,m1)
        real rx(m1)
        do 10 n=1,m1
        rx(n)=0
10      continue
        end
        
        subroutine zerm(mx,m1,m2)
        integer mx(m1,m2)
        do j=1,m2
        do i=1,m1
        mx(i,j)=0
	enddo
	enddo
        end
        
        subroutine zermr(rx,m1,m2)
        real rx(m1,m2)
        do j=1,m2
        do i=1,m1
        rx(i,j)=0
	enddo
	enddo
        end
	
	subroutine zermd(rx,m1,m2)
        double precision rx(m1,m2)
        do j=1,m2
        do i=1,m1
        rx(i,j)=0
	enddo
	enddo
        end
        
        subroutine eqm(j,k,m,n)
        integer j(m,n),k(m,n)
        do n1=1,n
        do m1=1,m
        j(m1,n1)=k(m1,n1)
	end do
        end do
        end
        
        
        real function rrand()
        double precision dseed,u
        save dseed
        data dseed /17395/
        call lrnd(dseed,u)
        rrand=real(u)
        end
        
        
        subroutine lrnd(dseed,u)
        double precision dseed, d31m1,u
        data d31m1 /2147483647/
        dseed=dmod(16087*dseed,d31m1)
        u=dseed/d31m1
        return
        end
	
	function rnorm(jd)
	rnorm=sqrt(-2*log(rrand()))*cos(6.283185*rrand())
	end 
	
        
        
        
        subroutine perm1(np,ns,tp)
        real tp(ns)
        j=np
11      rnd = rrand()
        k=int(j*rnd)+1
        tx=tp(j)
        tp(j)=tp(k)
        tp(k)=tx
        j=j-1
        if(j.gt.1) go to 11
        end
	