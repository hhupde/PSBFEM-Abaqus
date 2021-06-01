c ======================================================================
c User Subroutine UEL for sbfem 2D Elasticity
c All rights of reproduction or distribution in any form are reserved.
c By Yang Yang (PhD), yanghhu@foxmail.com
c ======================================================================
      include "EleCoeff2NodeEle.for"
      include "utilities.for"
c
      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1 PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2 KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     3 LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1 SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL),
     2 DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(*),
     3 JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4 PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)
      
      double precision ksc(mcrd,1)
      double precision RK_BD(ndofel,ndofel),EE0(ndofel,ndofel),EE1(ndofel,ndofel)
      double precision EE2(ndofel,ndofel)
      double precision E0(4,4),E1(4,4),E2(4,4),coord(2,2)
      double precision EE0H(ndofel,ndofel),EE1H(ndofel,ndofel),EE2H(ndofel,ndofel)
      double precision EE0H_inv(ndofel,ndofel),EE1H_T(ndofel,ndofel)
      double precision P(ndofel,ndofel),P_inv(ndofel,ndofel)
      double precision EE0P(ndofel,ndofel),EE1P(ndofel,ndofel),EE2P(ndofel,ndofel)
      double precision PP(2*ndofel,2*ndofel)
      double precision S1(ndofel,ndofel),S2(ndofel,ndofel)
      double precision S3(ndofel,ndofel),S4(ndofel,ndofel)
      double precision Zp(2*ndofel,2*ndofel),EEIT(ndofel,ndofel)
      double precision work(200),wr(2*ndofel),wi(2*ndofel),vl(2*ndofel,2*ndofel)
      double precision vr(2*ndofel,2*ndofel)
      double precision sortedArr(2*ndofel),eigvalue(ndofel),eigvector(2*ndofel,ndofel)
      double precision maxvalue(ndofel),eigvector_inv(ndofel,ndofel)
      double precision Kmatrix(ndofel,ndofel)
      double precision imatrix(ndofel,ndofel)
      double precision E,v
      double precision karea(nnode),kxc(nnode),kyc(nnode)
      integer ipiv(ndofel),id(2*ndofel)
      integer info
      integer jj,kk,k
c
      E=props(1)
      v=props(2)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Compute the coordinate of scaling centers | ksc
c     And translate the global coords to local coords
c
c     (1). Initial ksc
      if (kstep.eq.1.and.kinc.eq.1) then
      do i=1,mcrd
         ksc(i,1)=0.d0
      end do
      do i=1,nnode
         karea(i)=0.d0
      end do
      do i=1,nnode
         kxc(i)=0.d0
      end do
      do i=1,nnode
         kyc(i)=0.d0
      end do
c
c     (2). Compute ksc
      do i=1,nnode
        if (i.lt.nnode) then
         coord(1,1)=coords(1,i)
         coord(2,1)=coords(2,i)
         coord(1,2)=coords(1,i+1)
         coord(2,2)=coords(2,i+1)
         karea(i)=(coord(1,1)*coord(2,2)-coord(1,2)*coord(2,1))*0.5
         kxc(i)=(coord(1,1)+coord(1,2))/3.
         kyc(i)=(coord(2,1)+coord(2,2))/3.
        else
         coord(1,1)=coords(1,i)
         coord(2,1)=coords(2,i)
         coord(1,2)=coords(1,1)
         coord(2,2)=coords(2,1)
         karea(i)=(coord(1,1)*coord(2,2)-coord(1,2)*coord(2,1))*0.5
         kxc(i)=(coord(1,1)+coord(1,2))/3.
         kyc(i)=(coord(2,1)+coord(2,2))/3.
         end if
      end do
      do i=1,nnode
       ksc(1,1)=ksc(1,1)+karea(i)*kxc(i)/sum(karea)
       ksc(2,1)=ksc(2,1)+karea(i)*kyc(i)/sum(karea)
      end do
c
      do i=1,mcrd
         do j=1,nnode
            coords(i,j)=coords(i,j)-ksc(i,1)
         end do
      end do
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Compute EE0,EE1,EE2
c    (1). initial EE0,EE1,EE2
      do i=1,ndofel
         do j=1,ndofel
            EE0(i,j)=0.0D0
         end do
      end do
c     
      do i=1,ndofel
         do j=1,ndofel
            EE1(i,j)=0.0D0
         end do
      end do
c 
      do i=1,ndofel
         do j=1,ndofel
            EE2(i,j)=0.0D0
         end do
      end do
c    (2). calculating e0,e1,e2
       do i=1,nnode
        if (i.lt.nnode) then
         coord(1,1)=coords(1,i)
         coord(2,1)=coords(2,i)
         coord(1,2)=coords(1,i+1)
         coord(2,2)=coords(2,i+1)
c        
        call EleCoeff2NodeEle(coord,E,v,e0,e1,e2)
c
c     (3).   Assembly EE0,EE1,EE2
c
       do j=1,4
          do k=1,4
             jj=j+2*(i-1)
             kk=k+2*(i-1)
             EE0(jj,kk)=EE0(jj,kk)+e0(j,k)
          end do
        end do
c
       do j=1,4
          do k=1,4
             jj=j+2*(i-1)
             kk=k+2*(i-1)
             EE1(jj,kk)=EE1(jj,kk)+e1(j,k)
         end do
         end do
c
       do j=1,4
          do k=1,4
             jj=j+2*(i-1)
             kk=k+2*(i-1)
             EE2(jj,kk)=EE2(jj,kk)+e2(j,k)
         end do
         end do
c
      else
        coord(1,1)=coords(1,i)
        coord(2,1)=coords(2,i)
        coord(1,2)=coords(1,1)
        coord(2,2)=coords(2,1)
        call EleCoeff2NodeEle(coord,E,v,e0,e1,e2)
c
c        Assembly EE0,EE1,EE2
c
       do j=1,4
          do k=1,4
             jj=j+2*(i-1)
             kk=k+2*(i-1)
             if (kk.gt.ndofel) then
                 kk=kk-ndofel
             endif
             if (jj.gt.ndofel) then 
                 jj=jj-ndofel
             endif
             EE0(jj,kk)=EE0(jj,kk)+e0(j,k)
         end do
        end do
c
       do j=1,4
          do k=1,4
             jj=j+2*(i-1)
             kk=k+2*(i-1)
             if (kk.gt.ndofel) then
                 kk=kk-ndofel
             endif
             if (jj.gt.ndofel) then 
                 jj=jj-ndofel
             endif
             EE1(jj,kk)=EE1(jj,kk)+e1(j,k)
         end do
        end do
c
       do j=1,4
          do k=1,4
             jj=j+2*(i-1)
             kk=k+2*(i-1)
             if (kk.gt.ndofel) then
                 kk=kk-ndofel
             endif
             if (jj.gt.ndofel) then 
                 jj=jj-ndofel
             endif
             EE2(jj,kk)=EE2(jj,kk)+e2(j,k)
         end do
        end do
       endif
      end do
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Solver K matrix
c     Initiation matrix
      do i=1,ndofel
         do j=1,ndofel
            EE0H(i,j)=0.0d0
         end do
      end do
c
      do i=1,ndofel
         do j=1,ndofel
            EE1H(i,j)=0.0d0
         end do
      end do
c      
      do i=1,ndofel
         do j=1,ndofel
            EE2H(i,j)=0.0d0
         end do
      end do
c
      do i=1,ndofel
         do j=1,ndofel
            P(i,j)=0.0d0
         end do
      end do
c
      do i=1,ndofel
         P(i,i)=1.0d0/sqrt(abs(EE0(i,i)))
      end do
c
      EE0P=matmul(EE0,P)
      EE0H=matmul(P,EE0P)
      EE1P=matmul(EE1,P)
      EE1H=matmul(P,EE1P)
      EE2P=matmul(EE2,P)
      EE2H=matmul(P,EE2P)
c       
c     Initiation identity matrix
c
      do i=1,ndofel
         do j=1,ndofel
            imatrix(i,j)=0.0d0
         end do
      end do
      do i=1,ndofel
         imatrix(i,i)=1.0d0
      end do
c
c     calculate the inverse of EE0H
c
      call dgesv(ndofel,ndofel,EE0H,ndofel,ipiv,imatrix,ndofel,info)
      EE0H_inv=imatrix
c
c     calculate Zp
c
      EE1H_T=transpose(EE1H)
      S1=-matmul(EE0H_inv,EE1H_T)
      S2=EE0H_inv
      EEIT=matmul(EE0H_inv,EE1H_T)
      S3=EE2H-matmul(EE1H,EEIT)
      S4=matmul(EE1H,EE0H_inv)
c
      Zp(1:ndofel,1:ndofel)=S1
      Zp(1:ndofel,ndofel+1:2*ndofel)=S2
      Zp(ndofel+1:2*ndofel,1:ndofel)=S3
      Zp(ndofel+1:2*ndofel,ndofel+1:2*ndofel)=S4
c      
c     calculate eigenvalues and eignvectors
c
      wr=0.d0
      wi=0.d0
      vl=0.d0
      vr=0.d0
      call dgeev('N','V',2*ndofel,Zp,2*ndofel,wr,wi,vl,1,vr,
     1 2*ndofel,work,200,info)
c
      call sort(2*ndofel,wr,sortedArr,id,ndofel)
      eigvalue=sortedArr(1:ndofel)
      do i=1,2*ndofel
         do j=1,ndofel
            eigvector(i,j)=vr(i,id(j))
         end do
      end do
c
c     Initiation P_inv
c
      do i=1,ndofel
         do j=1,ndofel
            P_inv(i,j)=0.0d0
         end do
      end do
c
      do i=1,ndofel
         do j=1,ndofel
            P_inv(i,i)=1.0d0/P(i,i)
         end do
      end do
c
      PP(1:ndofel,1:ndofel)=P
      PP(1:ndofel,ndofel+1:2*ndofel)=0.0d0
      PP(ndofel+1:2*ndofel,1:ndofel)=0.0d0
      PP(ndofel+1:2*ndofel,ndofel+1:2*ndofel)=P_inv
c      
      eigvector=matmul(PP,eigvector)
c      
      eigvalue(ndofel-1:ndofel)=0.0d0
      eigvector(1:2*ndofel,ndofel-1:ndofel)=0.0d0
      do i=1,ndofel/2
         k=2*i-1
         eigvector(k,ndofel-1)=1.0d0
      end do
c
      do i=1,ndofel/2
         k=2*i
         eigvector(k,ndofel)=1.0d0
      end do
c
      do i=1,ndofel
         maxvalue(i)=maxval(abs(eigvector(1:ndofel-2,i)))
      end do
c
      do i=1,2*ndofel
         do j=1,ndofel
            eigvector(i,j)=eigvector(i,j)/maxvalue(j)
         end do
      end do
c
      do i=1,ndofel
         do j=1,ndofel
            imatrix(i,j)=0.0d0
         end do
      end do
      do i=1,ndofel
         imatrix(i,i)=1.0d0
      end do
c
c     calculate the inverse of eigvector
c
      call dgesv(ndofel,ndofel,eigvector(1:ndofel,1:ndofel),ndofel,ipiv,imatrix,ndofel,info)
c      print *, info
      eigvector_inv=imatrix
c
c     Solver K matrix
c
      Kmatrix=matmul(eigvector(ndofel+1:2*ndofel,1:ndofel),eigvector_inv)
      call STOREKmatrix(SVARS,NSVARS,Kmatrix,NDOFEL)
      else
      call READKmatrix(SVARS,NSVARS,Kmatrix,NDOFEL)
      end if

c     Initial rhs and amtrx
c
      do k1 = 1, NDOFEL
        do k2 = 1, NRHS
            rhs(k1,k2) = 0.d0
        end do
        do k2 = 1, NDOFEL
            amatrx(k2,k1) = 0.d0
        end do
      end do
c
c     Update K matrix
c
      do i = 1,ndofel
         do j = 1,ndofel
            amatrx(i,j) =  Kmatrix(i,j)
         end do
      end do
c
c     Update RHS
c
      do i = 1, ndofel
        do j = 1,ndofel
            rhs(i,1) = rhs(i,1) - amatrx(i,j) * U(j)
        end do
      end do 

      do i=1,ndofel
         do j=1,ndofel
            energy(2)=energy(2)+0.5d0*U(i)*Kmatrix(i,j)*U(j)
         end do
      end do
c
      return 
      end

c     Store stiff matrices
      SUBROUTINE STOREKmatrix(SVARS,NSVARS,Kmatrix,NDOFEL)
      
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION SVARS(NSVARS)
      DOUBLE PRECISION Kmatrix(NDOFEL,NDOFEL)
      
      do I = 1,NDOFEL
      SVARS((I-1)*NDOFEL+1:I*NDOFEL) = Kmatrix(I,1:NDOFEL)      
      end do
      
      RETURN
      
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      

c     Read stiff matrices      
      SUBROUTINE READKmatrix(SVARS,NSVARS,Kmatrix,NDOFEL)
      
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION SVARS(NSVARS)
      DOUBLE PRECISION Kmatrix(NDOFEL,NDOFEL)
c      
      do I = 1,NDOFEL
      Kmatrix(I,1:NDOFEL) = SVARS((I-1)*NDOFEL+1:I*NDOFEL)     
      end do

      RETURN
      
      END
