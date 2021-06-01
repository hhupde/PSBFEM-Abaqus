C ======================================================================
c
c User Subroutine EleCoeff2NodeEle for sbfem
c All rights of reproduction or distribution in any form are reserved.
c By Yang Yang (PhD), yanghhu@foxmail.com
c Parameter:
c          Coord: the coordinates of node | coord(dim,nnode)
c          E    : Young's modulus
c          v    : Possion's ratio
c Return:
c          the element coefficient matrices - e0, e1, e2
c
C ======================================================================
      subroutine EleCoeff2NodeEle(coord,E,v,E0,E1,E2)
      INCLUDE 'ABA_PARAM.INC'
c      implicit double precision(A-H,O-Z)
      double precision coord(2,2),C1(3,2),C2(3,2),D(3,3)
      double precision Q0(2,2),Q1(2,2),Q2(2,2),C1T(2,3),C2T(2,3)
      double precision Q0_temp(2,3),Q1_temp(2,3),Q2_temp(2,3)
      double precision E0(4,4),E1(4,4),E2(4,4)
      double precision meanX,meanY,deltaX,deltaY,area,Jb
      double precision E,v
c
c     Initiation Q0,Q1,Q2
c
      do i=1,2
         do j=1,2
            Q0(i,j)=0.0D0
         end do
      end do
c
      do i=1,2
         do j=1,2
            Q1(i,j)=0.0D0
         end do
      end do
c
      do i=1,2
         do j=1,2
            Q2(i,j)=0.0D0
         end do
      end do
c  
c     Initiation Q0_temp,Q1_temp,Q2_temp
      do i=1,2
         do j=1,3
            Q0_temp(i,j)=0.0D0
         end do
      end do
c
      do i=1,2
         do j=1,3
            Q1_temp(i,j)=0.0D0
         end do
      end do
c
      do i=1,2
         do j=1,3
            Q2_temp(i,j)=0.0D0
         end do
      end do
c
c     Initiation e0,e1,e2
      do i=1,4
         do j=1,4
            e0(i,j)=0.0D0
         end do
      end do
c   
      do i=1,4
         do j=1,4
            e1(i,j)=0.0D0
         end do
      end do
c
      do i=1,4
         do j=1,4
            e2(i,j)=0.0D0
         end do
      end do
c
      deltaX=coord(1,2)-coord(1,1)
      deltaY=coord(2,2)-coord(2,1)
      Jb=0.5d0*(coord(1,1)*coord(2,2)-coord(1,2)*coord(2,1))
      meanX=0.5d0*(coord(1,1)+coord(1,2))
      meanY=0.5d0*(coord(2,1)+coord(2,2))
      area=coord(1,1)*coord(2,2)-coord(1,2)*coord(2,1)
c
c      if (area.lt.1e-10) then 
c         write(*,*) "negative area"
c         write(*,*) "area",area
c         write(*,*) "coord",coord
c      end if
c
      C1(1,1)=deltaY*0.5d0
      C1(1,2)=0.0d0
      C1(2,1)=0.0d0
      C1(2,2)=-deltaX*0.5d0
      C1(3,1)=-deltaX*0.5d0
      C1(3,2)=deltaY*0.5d0
      C1T=transpose(C1)
c
      C2(1,1)=-meanY
      C2(1,2)=0.0d0
      C2(2,1)=0.0d0
      C2(2,2)=meanX
      C2(3,1)=meanX
      C2(3,2)=-meanY
      C2T=transpose(C2)
c
c     Constitutive matrix D
      D(1,1)=E/((1.0d0-v)*(1.0d0+v))
      D(1,2)=E*v/((1.0d0-v)*(1.0d0+v))
      D(1,3)=0.0d0
      D(2,1)=E*v/((1.0d0-v)*(1.0d0+v))
      D(2,2)=E/((1.0d0-v)*(1.0d0+v))
      D(2,3)=0.0d0
      D(3,1)=0.0d0
      D(3,2)=0.0d0
      D(3,3)=E/(2.0d0*(1.0d0+v))
c
c     Calculate Q0,Q1,Q2
      Q0_temp=0.5d0*(1.0d0/Jb)*matmul(C1T,D)
      Q0=matmul(Q0_temp,C1)
c
      Q1_temp=0.5d0*(1.0d0/Jb)*matmul(C2T,D)
      Q1=matmul(Q1_temp,C1)
c     
      Q2_temp=0.5d0*(1.0d0/Jb)*matmul(C2T,D)
      Q2=matmul(Q2_temp,C2)
c
c     Calculate E0,E1,E2
      e0(1:2,1:2)=2.0d0*Q0
      e0(1:2,3:4)=Q0
      e0(3:4,1:2)=Q0
      e0(3:4,3:4)=2.0d0*Q0
      e0=(2.0d0/3.0d0)*E0
c
      e1(1:2,1:2)=-1.0d0*Q1-(1.0d0/3.0d0)*Q0
      e1(1:2,3:4)=-1.0d0*Q1+(1.0d0/3.0d0)*Q0
      e1(3:4,1:2)=Q1+(1.0d0/3.0d0)*Q0
      e1(3:4,3:4)=Q1-(1.0d0/3.0d0)*Q0
c
      e2(1:2,1:2)=Q2+(1.0d0/3.0d0)*Q0
      e2(1:2,3:4)=-1.0d0*Q2-(1.0d0/3.0d0)*Q0
      e2(3:4,1:2)=-1.0d0*Q2-(1.0d0/3.0d0)*Q0
      e2(3:4,3:4)=Q2+(1.0d0/3.0d0)*Q0
c
      return
      end
