cccccccccc FORTRAN subroutine rlbin.f cccccccccc

c Obtains bin counts for univariate regression data
c via the linear binning strategy. If "trun=0" then
c weight from end observations is given to corresponding
c end grid points. If "trun=1" then end observations
c are truncated.

c Last changed: 06/02/95

      subroutine rlbin(X,Y,n,a,b,M,trun,xcounts,ycounts)     
      double precision X(*),Y(*),a,b,xcounts(*),ycounts(*),lxi,delta,rem
      integer n,M,i,li,trun

c     Initialize grid counts to zero

      do 10 i=1,M
         xcounts(i) = dble(0)
         ycounts(i) = dble(0)
10    continue

      delta = (b-a)/(M-1)
      do 20 i=1,n
         lxi = ((X(i)-a)/delta) + 1

c        Find integer part of "lxi"

         li = lxi 
         rem = lxi - li
         if (li.ge.1.and.li.lt.M) then
            xcounts(li) = xcounts(li) + (1-rem)
            xcounts(li+1) = xcounts(li+1) + rem
            ycounts(li) = ycounts(li) + (1-rem)*y(i)
            ycounts(li+1) = ycounts(li+1) + rem*y(i)
         elseif (li.lt.1.and.trun.eq.0) then
            xcounts(1) = xcounts(1) + 1
            ycounts(1) = ycounts(1) + y(i)
         elseif (li.ge.M) then
            if (li.eq.M.or.trun.eq.0) then 
               xcounts(M) = xcounts(M) + 1
               ycounts(M) = ycounts(M) + y(i)
            endif
         endif
20    continue

      return
      end

cccccccccc End of rlbin.f cccccccccc
