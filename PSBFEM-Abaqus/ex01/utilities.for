      subroutine sort(n,a,tempArr,id,ndofel)
      INCLUDE 'ABA_PARAM.INC' 
      double precision a(2*ndofel),tempArr(2*ndofel)
      double precision min
      integer id(2*ndofel)
      tempArr=a
      min=MinVal(a)-1.0d0
c     sort
      do i=1,n-1
         do j=i+1,n
           if (tempArr(i) .lt. tempArr(j)) then
             temp = tempArr(i)
             tempArr(i) = tempArr(j)
             tempArr(j) = temp
           end if
         end do
      end do
c    obtain index
      do i =1,n
         do j=1,n
          if(tempArr(i).eq.a(j)) then
             id(i)=j
             a(j)=min
             exit
          end if
         end do
      end do
      return
      end
